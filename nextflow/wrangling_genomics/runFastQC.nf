untrimmed_FastQs = Channel.fromPath("$HOME/working/wrangling-genomics/data/untrimmed_fastq/*.fastq.gz")

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