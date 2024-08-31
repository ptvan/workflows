!#/bin/sh 

export BT2IDX=/Users/ptv/working/Databases

## download data
# NOTE: below are 1/400 subsampled versions
# full versions are in SRA (ENCSR356KRQ)

wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF341MYG.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF106QGY.subsampled.400.fastq.gz

wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF248EJF.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF368TYI.subsampled.400.fastq.gz

wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF641SFZ.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF751XTV.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF927LSG.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF859BDM.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF193RRC.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF366DFI.subsampled.400.fastq.gz

wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF031ARQ.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF590SYZ.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF734PEQ.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF007USV.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF886FSC.subsampled.400.fastq.gz
wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF573UXK.subsampled.400.fastq.gz


## download databases
wget https://storage.googleapis.com/encode-pipeline-genome-data/hg38/hg38.blacklist.bed.gz
wget https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/hg38_caper.tsv
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

## build bowtie2 index
# NOTE: this took ~2 hrs, maybe just download indices from hg38_caper.tsv location ?
bowtie2-build GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz GCRh38_ATACseq

## trim sequencing adaptors
fastp -i ENCFF641SFZ.subsampled.400.fastq.gz \
      -I ENCFF031ARQ.subsampled.400.fastq.gz \
      -o sample1_R1.trimmed.fastq.gz \
      -O sample1_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample1.fastp.json \
      -h sample1.fastp.html

## perform alignment against GRCh38
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample1_R1.trimmed.fastq.gz -2 sample1_R2.trimmed.fastq.gz | samtools view -bS - > sample1.bam

## sort and index BAMs
samtools sort sample1.bam -o sample1_sorted.bam
samtools index sample1_sorted.bam -o sample1_sorted.bai

## generate mapping stats
samtools idxstats sample1_sorted.bam > sample1_sorted.idxstats
grep "chrM" sample1_sorted.idxstats
samtools flagstat sample1_sorted.bam > sample1_sorted.flagstat

## remove mitochondrial reads
samtools view -h sample1_sorted.bam | grep -v chrM | samtools sort -O bam -o sample1_sorted.rmChrM.bam -T .
samtools index sample1_sorted.rmChrM.bam -o sample1_sorted.rmChrM.bai

## check if there is a readgroup (@RG tag)
samtools view -H sample1_sorted.rmChrM.bam | grep '@RG'

## ... if not, add @RG tag to BAM
samtools addreplacerg -r "@RG\tID:RG1\tSM:Sample1\tPL:Illumina\tLB:Library.fa" -o sample1_sorted.rmChrM.RGadded.bam sample1_sorted.rmChrM.bam

## markDuplicates
java -jar ~/working/packages/picard.jar MarkDuplicates QUIET=true INPUT=sample1_sorted.rmChrM.RGadded.bam OUTPUT=sample1_sorted.rmChrM.RGadded.dupesmarked.bam METRICS_FILE=sample1.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

## remove duplicates
samtools view -q 30 -c sample1_sorted.rmChrM.RGadded.dupesmarked.bam

samtools view -h -b -f 2 -F 1548 -q 30 sample1_sorted.rmChrM.RGadded.dupesmarked.bam | samtools sort -o sample1_sorted.rmChrM.RGadded.dupesmarked.filtered.bam
samtools index sample1_sorted.rmChrM.RGadded.dupesmarked.filtered.bam -o sample1_sorted.rmChrM.RGadded.dupesmarked.filtered.bai

## remove reads within hg38 blacklist regions
bedtools intersect -nonamecheck -v -abam sample1_sorted.rmChrM.RGadded.dupesmarked.filtered.bam -b hg38.blacklist.bed.gz > sample1.tmp.bam

## sort and index
samtools sort -O bam -o sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.bam sample1.tmp.bam
samtools index sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.bam -o sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.bai
rm sample1.tmp.bam

## shift read coordinates (optional)
alignmentSieve \
--verbose \
--ATACshift \
--blackListFileName hg38.blacklist.bed.gz \
--bam sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.bam \
-o sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bam

samtools sort -O bam -o sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bam sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bam
samtools index sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bam -o sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bai


## bam to bigwig
# 2862010578 is the effective genome size for GRCh38 when using 150bp reads and including only regions which are uniquely mappable

bamCoverage \
--numberOfProcessors 8 \
--binSize 10 \
--normalizeUsing BPM \
--effectiveGenomeSize 2862010578 \
--bam sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bam \
-o sample1_coverage_BPM.bw

## bam to BEDPE
macs3 randsample -i sample1_sorted.rmChrM.RGadded.dupesmarked.blacklist-filtered.shifted.bam -f BAMPE -p 100 -o sample1.bed

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
-t sample1.bed \
--outdir macs3/sample1 2> macs3.log