!#/bin/sh 

export BT2IDX=/Users/ptv/working/Databases

## perform alignment against GRCh38
bowtie2 \
--local \
--very-sensitive \
--no-mixed \
--no-discordant \
-I 25 \
-X 700 \
-x $BT2IDX/GCRh38_ATACseq \
-1 sample1_R1.trimmed.fastq.gz \
-2 sample1_R2.trimmed.fastq.gz | samtools view -bS - > sample1.bam

## sort and index BAMs
samtools sort sample1.bam -o sample1.tmp.bam
samtools index sample1.tmp.bam -o sample1.bai
mv sample1.tmp.bam sample1.bam

## generate mapping stats
samtools idxstats sample1.bam > sample1.idxstats
grep "chrM" sample1.idxstats
samtools flagstat sample1.bam > sample1.flagstat

## remove mitochondrial reads
samtools view -h sample1.bam | grep -v chrM | samtools sort -O bam -o sample1.noChrM.bam -T .
samtools index sample1.noChrM.bam -o sample1.noChrM.bai

## check if there is a readgroup (@RG tag)
samtools view -H sample1.noChrM.bam | grep '@RG'

## ... if not, add @RG tag to BAM
samtools addreplacerg -r "@RG\tID:RG1\tSM:Sample1\tPL:Illumina\tLB:Library.fa" -o sample1.noChrM.RGadded.bam sample1.noChrM.bam

## markDuplicates
# java -jar ~/working/packages/picard.jar MarkDuplicates QUIET=true INPUT=sample1.noChrM.RGadded.bam OUTPUT=sample1.noChrM.RGadded.noDuplicates.bam METRICS_FILE=sample1.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

## remove duplicates
#samtools view -q 30 -c sample1.noChrM.RGadded.noDuplicates.bam

samtools view -h -b -f 2 -F 1548 -q 30 sample1.noChrM.RGadded.bam | samtools sort -o sample1.noChrM.RGadded.noDuplicates.bam
samtools index sample1.noChrM.RGadded.noDuplicates.bam -o sample1.noChrM.RGadded.noDuplicates.bai

## remove reads within hg38 blacklist regions
bedtools intersect -nonamecheck -v -abam sample1.noChrM.RGadded.noDuplicates.bam -b ../hg38.blacklist.bed.gz > sample1.tmp.bam

samtools sort -O bam -o sample1.noChrM.RGadded.noDuplicates.blacklisted.bam sample1.tmp.bam
samtools index sample1.noChrM.RGadded.noDuplicates.blacklisted.bam -o sample1.noChrM.RGadded.noDuplicates.blacklisted.bai
rm sample1.tmp.bam

## shift read coordinates (optional)
alignmentSieve \
--verbose \
--ATACshift \
--blackListFileName ../hg38.blacklist.bed.gz \
--bam sample1.noChrM.RGadded.noDuplicates.blacklisted.bam \
-o sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam

samtools sort -O bam -o sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam
samtools index sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam -o sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bai

## bam to bigwig
# 2862010578 is the effective genome size for GRCh38 when using 150bp reads and including only regions which are uniquely mappable
bamCoverage \
--numberOfProcessors 8 \
--binSize 10 \
--normalizeUsing BPM \
--effectiveGenomeSize 2862010578 \
--bam sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam \
-o sample1_coverage_BPM.bw

## bam to BEDPE
# macs3 randsample -i sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam -f BAMPE -p 100 -o sample1.bed
bedtools bamtobed -i sample1.noChrM.RGadded.noDuplicates.blacklisted.shifted.bam -bedpe > sample1.bedPE

## call peaks
macs3 callpeak \
-f BEDPE \
--nomodel \
--shift -37 \
--extsize 73 \
-g 2862010578 \
-B --broad \
--keep-dup all \
--cutoff-analysis -n sample1 \
-t sample1.bedPE \
--outdir ./ 2> macs3.log
