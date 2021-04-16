nextflow.enable.dsl=2

params.suffix = "_{1,2}.trim.fastq"
params.alignmentProgram = "/usr/bin/bwa"
params.alignmentMode = "mem"
params.referenceGenome = "$HOME/working/wrangling-genomics/data/ref_genome/ecoli_rel606.fasta"
params.help = false

if (params.help || params.input == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run align.nf --input <> --output <>
    
    Required Arguments:
      --input        Full path to .fasta.gz file to trim
      --output       Folder to place trimmed files
    
    """.stripIndent()
}


Channel.fromFilePairs("${params.input}**${params.suffix}").set{in_fastq}

process align {
    publishDir "${params.output}", mode:"copy", overwrite: true

    input:
    tuple val(name), file(in_fastq)

    output:
    file("*.sam") 

    script:
    """
    ${params.alignmentProgram} ${params.alignmentMode} ${params.referenceGenome} \
              ${in_fastq.get(0)} ${in_fastq.get(1)} > ${in_fastq.get(0).getBaseName().take(in_fastq.get(0).name.lastIndexOf('_1'))}.sam
    """
}

workflow {

    main:
        align(in_fastq)
    
    emit:
        align.out    

}