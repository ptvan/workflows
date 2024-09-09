#!/bin/bash

## this ATACseq workflow draws from several sources:
## https://github.com/ENCODE-DCC/atac-seq-pipeline
## https://github.com/CebolaLab/ATAC-seq
## https://vallierlab.wixsite.com/pipelines/atac-seq

## it takes a directory of BAMs and outputs all output and intermediate files in another directory
## example data for a short (~15 minutes) run is preceded by ATACseq_download_and_setup.sh
## which fetches public data, creates FASTQs and performs alignment prior to this workflow

## PREREQUISITES:
## - bowtie2
## - samtools
## - bedtools
## - picard
## - deeptools
## - MACS3

## after peak-calling, additional analyses can be performed outside this workflow
## eg. R/ATACSeq_postpeakcalling_workflow.R in this repo

while getopts ":i:o:" opt; do
      case $opt in
      i) arg_1="$OPTARG"
      ;;
      o) arg_2="$OPTARG"
      ;;
      esac
done

if [ -z $arg_1 ]
then 
      printf " -i cannot be empty !!! Exiting... \n "
      printf "\n USAGE: $0 -i <directory with BAM file(s)> -o <output directory> \n" 
      exit
else
      printf "checking input BAM directory ${arg_1} ....  OK !\n"
fi

if [ -z $arg_2 ]
then 
      printf " -o cannot be empty !!! Exiting... \n "
      printf "\n USAGE: $0 -i <directory with BAM file(s)> -o <output directory> \n" 
      exit
else
      outputdir=$arg_2
      printf "checking output directory ${outputdir} ... OK !\n"
      
fi

printf "\n RUNNING ATACseq WORKFLOW... \n\n"

for inputpath in `ls $arg_1/*.bam`
do 
      filename="${inputpath##*/}"
      samplename="${filename%.bam}"
      samplenumber="${samplename:6:1}"
      printf "processing $samplename in $inputpath: \n"
      
      printf " >>> index and sort ... \n"
      samtools sort $inputpath -o $outputdir/$samplename.sorted.bam
      samtools index $outputdir/$samplename.sorted.bam -o $outputdir/$samplename.sorted.bai

      printf " >>> generate mapping stats ... \n"
      samtools idxstats $outputdir/$samplename.sorted.bam > $outputdir/$samplename.sorted.idxstats
      samtools flagstats $outputdir/$samplename.sorted.bam > $outputdir/$samplename.sorted.flagstats
      
      printf " >>> removing mitochondrial reads \n"
      samtools view -h $outputdir/$samplename.sorted.bam | grep -v chrM | samtools sort -O bam -o $outputdir/$samplename.sorted.rmChrM.bam -T .
      samtools index $outputdir/$samplename.sorted.rmChrM.bam -o $outputdir/$samplename.sorted.rmChrM.bai

      printf " >>> checking for readgroup ..."
      if samtools view -H $outputdir/$samplename.sorted.rmChrM.bam | grep '@RG'
      then
            printf "found @RG tag \n"
      else 
            printf "NO READGROUP FOUND !!! "
            samtools addreplacerg -r "@RG\tID:RG$samplenumber\tSM:$samplename\tPL:Illumina\tLB:Library.fa" -o $outputdir/$samplename.sorted.rmChrM.RGadded.bam $outputdir/$samplename.sorted.rmChrM.bam
            mv $outputdir/$samplename.sorted.rmChrM.RGadded.bam $outputdir/$samplename.sorted.rmChrM.bam
            printf "Added @RG tag to BAM \n"  
      fi

      printf " >>> marking & removing duplicate reads ..."
      java -jar ~/working/packages/picard.jar MarkDuplicates \
            QUIET=true \
            INPUT=$outputdir/$samplename.sorted.rmChrM.bam \
            OUTPUT=$outputdir/$samplename.sorted.rmChrM.nodupes.bam \
            METRICS_FILE=$outputdir/$samplename.dup.metrics \
            REMOVE_DUPLICATES=true \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT TMP_DIR=. > $outputdir/$samplename.MarkDuplicates.log 2> /dev/null 
      printf "\n"      

      printf " >>> removing blacklisted regions ... \n"
      bedtools intersect -nonamecheck -v -abam $outputdir/$samplename.sorted.rmChrM.nodupes.bam -b hg38.blacklist.bed.gz > $outputdir/$samplename.tmp.bam
      samtools sort -O bam -o $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.bam $outputdir/$samplename.tmp.bam
      samtools index $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.bam -o $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.bai
      rm $outputdir/$samplename.tmp.bam

      printf " >>> shifting coordinates ... \n"
      alignmentSieve \
            --verbose \
            --ATACshift \
            --blackListFileName hg38.blacklist.bed.gz \
            --bam $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.bam \
            -o $outputdir/$samplename.tmp.bam > $outputdir/$samplename.alignmentSieve.log 2> /dev/null 
      samtools sort -O bam -o $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.shifted.bam $outputdir/$samplename.tmp.bam
      samtools index $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.shifted.bam -o $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.shifted.bai
      rm $outputdir/$samplename.tmp.bam

      printf " >>> generating coverage BigWig ... \n"
      # 2862010578 is the effective genome size for GRCh38 when using 150bp reads and including only regions which are uniquely mappable
      bamCoverage \
            --numberOfProcessors 8 \
            --binSize 10 \
            --normalizeUsing BPM \
            --effectiveGenomeSize 2862010578 \
            --bam $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.shifted.bam \
            -o $outputdir/$samplename.coverageBPM.bw > $outputdir/$samplename.bamCoverage.log 2> /dev/null 

      printf " >>> converting BAM to BEDPE ... \n"
      macs3 randsample \
            -i $outputdir/$samplename.sorted.rmChrM.nodupes.filtered.shifted.bam \
            -f BAMPE \
            -p 100 \
            -o $outputdir/$samplenamePE.bed > $outputdir/$samplename.randsample.log 2> /dev/null 

      printf " >>> calling peaks ... \n"
      macs3 callpeak \
            -f BEDPE \
            --nomodel \
            --shift -37 \
            --extsize 73 \
            -g 2862010578 \
            -B --broad \
            --keep-dup all \
            --cutoff-analysis -n $samplename \
            -t $outputdir/$samplenamePE.bed \
            --outdir ./$outputdir 2> $outputdir/$samplename.callpeak.log

      printf "$samplename COMPLETED SUCCESSFULLY ! \n\n"
done