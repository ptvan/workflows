nextflow.enable.dsl=2

params.suffix = "_{1,2}.fastq.gz"
params.trimProgram = "$HOME/working/trimmomatic/trimmomatic-0.39.jar"
params.trimMode = "PE"

adapterFile = "$HOME/working/trimmomatic/adapters/NexteraPE-PE.fa"

Channel.fromFilePairs("${params.input}**${params.suffix}").set{in_fastq}

process trimAdaptors {
//    publishDir "$HOME/working/wrangling-genomics/data/trimmed_fastq_small/", mode:"copy", overwrite: true

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
              SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${adapterFile}:2:40:15 
    """
}

workflow {

    main:
        trimAdaptors(in_fastq)
    
    emit:
        trimAdaptors.out    

}