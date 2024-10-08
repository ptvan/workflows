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

## trim sequencing adaptors & rename to R1/R2 convention
fastp -i ENCFF641SFZ.subsampled.400.fastq.gz \
      -I ENCFF031ARQ.subsampled.400.fastq.gz \
      -o sample1_R1.trimmed.fastq.gz \
      -O sample1_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample1.fastp.json \
      -h sample1.fastp.html

fastp -i ENCFF751XTV.subsampled.400.fastq.gz \
      -I ENCFF590SYZ.subsampled.400.fastq.gz \
      -o sample2_R1.trimmed.fastq.gz \
      -O sample2_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample2.fastp.json \
      -h sample2.fastp.html

fastp -i ENCFF927LSG.subsampled.400.fastq.gz \
      -I ENCFF734PEQ.subsampled.400.fastq.gz \
      -o sample3_R1.trimmed.fastq.gz \
      -O sample3_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample3.fastp.json \
      -h sample3.fastp.html

fastp -i ENCFF859BDM.subsampled.400.fastq.gz \
      -I ENCFF007USV.subsampled.400.fastq.gz \
      -o sample4_R1.trimmed.fastq.gz \
      -O sample4_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample4.fastp.json \
      -h sample4.fastp.html

fastp -i ENCFF193RRC.subsampled.400.fastq.gz \
      -I ENCFF886FSC.subsampled.400.fastq.gz \
      -o sample5_R1.trimmed.fastq.gz \
      -O sample5_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample5.fastp.json \
      -h sample5.fastp.html

fastp -i ENCFF366DFI.subsampled.400.fastq.gz \
      -I ENCFF573UXK.subsampled.400.fastq.gz \
      -o sample6_R1.trimmed.fastq.gz \
      -O sample6_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample6.fastp.json \
      -h sample6.fastp.html

fastp -i ENCFF341MYG.subsampled.400.fastq.gz \
      -I ENCFF248EJF.subsampled.400.fastq.gz \
      -o sample7_R1.trimmed.fastq.gz \
      -O sample7_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample7.fastp.json \
      -h sample7.fastp.html

fastp -i ENCFF106QGY.subsampled.400.fastq.gz \
      -I ENCFF368TYI.subsampled.400.fastq.gz \
      -o sample8_R1.trimmed.fastq.gz \
      -O sample8_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      -j sample8.fastp.json \
      -h sample8.fastp.html

bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample1_R1.trimmed.fastq.gz -2 sample1_R2.trimmed.fastq.gz | samtools view -bS - > input/sample1.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample2_R1.trimmed.fastq.gz -2 sample2_R2.trimmed.fastq.gz | samtools view -bS - > input/sample2.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample3_R1.trimmed.fastq.gz -2 sample3_R2.trimmed.fastq.gz | samtools view -bS - > input/sample3.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample4_R1.trimmed.fastq.gz -2 sample4_R2.trimmed.fastq.gz | samtools view -bS - > input/sample4.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample5_R1.trimmed.fastq.gz -2 sample5_R2.trimmed.fastq.gz | samtools view -bS - > input/sample5.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample6_R1.trimmed.fastq.gz -2 sample6_R2.trimmed.fastq.gz | samtools view -bS - > input/sample6.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample7_R1.trimmed.fastq.gz -2 sample7_R2.trimmed.fastq.gz | samtools view -bS - > input/sample7.bam
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $BT2IDX/GCRh38_ATACseq -1 sample8_R1.trimmed.fastq.gz -2 sample8_R2.trimmed.fastq.gz | samtools view -bS - > input/sample8.bam