nextflow.enable.dsl=2

params.suffix = ".fasta.*"
params.help = false

if (params.help || params.file == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run indexReferenceGenome.nf --file <> --output <>
    
    Required Arguments:
      --input        Full path to .fasta file to index
      --output       Folder to place index files
    
    """.stripIndent()
}

process indexReferenceGenome {
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    file(ref)

    output:
    file "*${params.suffix}"

    script:
    """
    bwa index ${ref}
    """
}

workflow {
    input_ch = Channel.fromPath(
        "${params.file}"
    )

    indexReferenceGenome(
        input_ch
    )

}