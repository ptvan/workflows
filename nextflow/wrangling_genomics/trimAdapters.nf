nextflow.enable.dsl=2

params.suffix = "_{1,2}.fastq.gz"
params.trimProgram = "$HOME/working/trimmomatic/trimmomatic-0.39.jar"
params.trimMode = "PE"
params.adapterFile = "$HOME/working/trimmomatic/adapters/NexteraPE-PE.fa"
params.help = false

if (params.help || params.input == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run trimAdapters.nf --input <> --output <>
    
    Required Arguments:
      --input        Full path to .fasta.gz file to trim
      --output       Folder to place trimmed files
    
    """.stripIndent()
}


Channel.fromFilePairs("${params.input}**${params.suffix}").set{in_fastq}

process trimAdaptors {
    publishDir "${params.output}", mode:"copy", overwrite: true

    input:
    tuple val(name), file(in_fastq)

    output:
    file("*.trim.fastq") 

    script:
    """
    java -jar ${params.trimProgram} ${params.trimMode} -threads 4 \
               ${in_fastq.get(0)} \
              ${in_fastq.get(1)}  \
              ${in_fastq.get(0).simpleName}.trim.fastq ${in_fastq.get(0).simpleName}un.trim.fastq \
              ${in_fastq.get(1).simpleName}.trim.fastq ${in_fastq.get(1).simpleName}un.trim.fastq \
              SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${params.adapterFile}:2:40:15 
    """
}

workflow {

    main:
        trimAdaptors(in_fastq)
    
    emit:
        trimAdaptors.out    

}