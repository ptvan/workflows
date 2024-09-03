#!/bin/bash

for filename in `ls sample[1-9].bam`
do 
      samplename="${filename%.bam}"
      samplenumber="${samplename:6:1}"
      printf "processing $samplename: \n"
      
      printf " >>> index and sort \n"
      samtools sort $filename -o $samplename.sorted.bam
      samtools index $samplename.sorted.bam -o $samplename.sorted.bai

      printf " >>> generate mapping stats \n"
      samtools idxstats $samplename.sorted.bam > $samplename.sorted.idxstats
      samtools flagstats $samplename.sorted.bam > $samplename.sorted.flagstats
      
      printf " >>> removing mitochondiral reads \n"
      samtools view -h $samplename.sorted.bam | grep -v chrM | samtools sort -O bam -o $samplename.sorted.rmChrM.bam -T .
      samtools index $samplename.sorted.rmChrM.bam -o $samplename.sorted.rmChrM.bai

      printf " >>> checking for readgroup ..."
      if samtools view -H sample1.bam | grep '@RG'
      then
            printf "found @RG tag \n"
      else 
            printf "NO READGROUP FOUND !!! "
            samtools addreplacerg -r "@RG\tID:RG$samplenumber\tSM:$samplename\tPL:Illumina\tLB:Library.fa" -o $samplename.sorted.rmChrM.RGadded.bam $samplename.sorted.rmChrM.bam
            mv $samplename.sorted.rmChrM.RGadded.bam $samplename.sorted.rmChrM.bam
            printf "Added @RG tag \n"  
      fi

      printf " >>> marking & removing duplicate reads ..."
      java -jar ~/working/packages/picard.jar MarkDuplicates \
            QUIET=true \
            INPUT=$samplename.sorted.rmChrM.bam \
            OUTPUT=$samplename.sorted.rmChrM.nodupes.bam \
            METRICS_FILE=$samplename.dup.metrics \
            REMOVE_DUPLICATES=true \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT TMP_DIR=. > $samplename.MarkDuplicates.log 2> /dev/null 
      printf "\n"      

      printf " >>> removing blacklisted regions \n"
      bedtools intersect -nonamecheck -v -abam $samplename.sorted.rmChrM.nodupes.bam -b hg38.blacklist.bed.gz > $samplename.tmp.bam
      samtools sort -O bam -o $samplename.sorted.rmChrM.nodupes.filtered.bam $samplename.tmp.bam
      samtools index $samplename.sorted.rmChrM.nodupes.filtered.bam -o $samplename.sorted.rmChrM.nodupes.filtered.bai
      rm $samplename.tmp.bam

      printf " >>> shifting coordinates \n"

      alignmentSieve \
            --verbose \
            --ATACshift \
            --blackListFileName hg38.blacklist.bed.gz \
            --bam $samplename.sorted.rmChrM.nodupes.filtered.bam \
            -o $samplename.tmp.bam > $samplename.alignmentSieve.log 2> /dev/null 

      samtools sort -O bam -o $samplename.sorted.rmChrM.nodupes.filtered.shifted.bam $samplename.tmp.bam
      samtools index $samplename.sorted.rmChrM.nodupes.filtered.shifted.bam -o $samplename.sorted.rmChrM.nodupes.filtered.shifted.bai

      printf " >>> generating coverage BigWig \n"
      # 2862010578 is the effective genome size for GRCh38 when using 150bp reads and including only regions which are uniquely mappable
      bamCoverage \
            --numberOfProcessors 8 \
            --binSize 10 \
            --normalizeUsing BPM \
            --effectiveGenomeSize 2862010578 \
            --bam $samplename.sorted.rmChrM.nodupes.filtered.shifted.bam \
            -o $samplename.coverageBPM.bw > $samplename.bamCoverage.log 2> /dev/null 

      printf " >>> converting BAM to BEDPE \n"
      macs3 randsample \
            -i $samplename.sorted.rmChrM.nodupes.filtered.shifted.bam \
            -f BAMPE \
            -p 100 \
            -o $samplenamePE.bed > $samplename.randsample.log 2> /dev/null 

      printf " >>> calling peaks \n"
      macs3 callpeak \
            -f BEDPE \
            --nomodel \
            --shift -37 \
            --extsize 73 \
            -g 2862010578 \
            -B --broad \
            --keep-dup all \
            --cutoff-analysis -n $samplename \
            -t $samplenamePE.bed \
            --outdir ./ 2> $samplename.callpeak.log

      printf "$samplename COMPLETED SUCCESSFULLY ! \n\n"
done