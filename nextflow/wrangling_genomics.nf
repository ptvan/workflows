#!/usr/bin/env nextflow
 
params.input = "$HOME/working/wrangling-genomics/data/untrimmed_fastq/*.fastq.gz"

untrimmed_FastQs = Channel.fromPath(params.input)

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
