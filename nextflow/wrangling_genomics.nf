#!/usr/bin/env nextflow
// based on https://datacarpentry.org/wrangling-genomics/

params.input = "$HOME/working/wrangling-genomics/data/untrimmed_fastq/*.fastq.gz"
untrimmed_FastQs = Channel.fromPath(params.input)

read_pairs = Channel.fromFilePairs("$HOME/working/wrangling-genomics/data/untrimmed_fastq/*_[1,2].fastq.gz", type: 'file')

process runFastQC {
    publishDir "$HOME/working/wrangling-genomics/data/untrimmed_fastq/", mode:"copy", overwrite: true

    input:
    file untrimmed_fastq from untrimmed_FastQs

    output:
    file("*.html") into fastQCs

    script:
    """
    fastqc ${untrimmed_fastq}
    """

}

process trimAdaptors {
    publishDir "$HOME/working/wrangling-genomics/data/trimmed_fastq_small/", mode:"copy", overwrite: true

    input:
    set sample, file(in_fastq) from read_pairs

    output:
    file("*.trim.fastq.gz") into trimmed_FastQs

    script:
    """
    java -jar $HOME/working/trimmomatic/trimmomatic-0.39.jar PE -threads 4 ${in_fastq.get(0)}.fastq ${in_fastq.get(1)}.fastq  \
              ${in_fastq.get(0)}.trimmed.fastq ${in_fastq.get(0)}un.trimmed.fastq \
              ${in_fastq.get(1)}.trimmed.fastq ${in_fastq.get(1)}un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
    """
}
