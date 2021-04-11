#!/usr/bin/env nextflow
// based on https://datacarpentry.org/wrangling-genomics/

params.input = "$HOME/working/wrangling-genomics/data/untrimmed_fastq/*.fastq.gz"
untrimmed_FastQs = Channel.fromPath(params.input)

read_pairs = Channel.fromFilePairs("$HOME/working/wrangling-genomics/data/untrimmed_fastq/*_[1,2].fastq.gz", type: 'file')

adapterFile = "$HOME/working/trimmomatic/adapters/NexteraPE-PE.fa"

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
    file("*.trim.fastq") into trimmed_FastQs

    script:
    """
    java -jar $HOME/working/trimmomatic/trimmomatic-0.39.jar PE -threads 4 \
              ${in_fastq.get(0)} \
              ${in_fastq.get(1)}  \
              ${in_fastq.get(0).simpleName}.trim.fastq ${in_fastq.get(0).simpleName}un.trim.fastq \
              ${in_fastq.get(1).simpleName}.trim.fastq ${in_fastq.get(1).simpleName}un.trim.fastq \
              SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${adapterFile}:2:40:15 
    """
}
