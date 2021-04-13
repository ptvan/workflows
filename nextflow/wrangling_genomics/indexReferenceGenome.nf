referenceGenomeFile = "$HOME/working/wrangling-genomics/data/ref_genome/ecoli_rel606.fasta"

referenceGenomeIndex = "$HOME/working/wrangling-genomics/data/ref_genome/ecoli_rel606.fasta.bwt"

process indexReferenceGenome {
    publishDir "$HOME/working/wrangling-genomics/data/ref_genome/", mode:"copy", overwrite: true

    input:
    file x from referenceGenomeFile

    output:
    file("${x}") into genomeIndex

    script:
    """
    bwa index ${referenceGenomeFile}
    """
}

