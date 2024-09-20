nextflow.enable.dsl=2

params.samSuffix = ".sam"
params.samtoolsProgram = "/usr/local/bin/samtools"
params.samtoolsConvertParams = "view -S -b "
params.samtoolsSortParams = "sort -o"
params.bamSuffix = ".bam"
params.help = false

if (params.help || params.input == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run samToSortedBam.nf --input <> --output <>
    
    Required Arguments:
      --input        Full path to .sam file to convert and index to .bam
      --output       Folder to place convert indexed files
    
    NOTE: this pipeline requires GNU parallel 
    """.stripIndent()
}

Channel.fromFile("${params.input}**${params.samSuffix}").set{in_sam}
Channel.fromFile("${params.input}**${params.bamSuffix}").set{in_bam}


process convertToBam {
    publishDir "${params.output}", mode:"copy", overwrite: true

    input:
      tuple val(name), file(in_sam)
   
    output:
      path '*.bam', emit: unsorted_bam_ch

    script:
    """
    ${params.samtoolsProgram} ${params.samtoolsConvertParams} \
              ${in_sam.get(0)} > ${in_sam.get(0).getBaseName()}.bam
    ls -lahtr
    """
}

process sortBam {
    publishDir "${params.output}", mode:"copy", overwrite: true

    input:
      tuple val(name), file(in_bam)
   
    output:
      path '*.bam', emit: sorted_bam_ch

    script:
    """
    ${params.samtoolsProgram} ${params.samtoolsSortParams} \
              ${in_bam.get(0)} > ${in_bam.get(0).getBaseName()}.sorted.bam
    ls -lahtr
    """
}

workflow {
    convertToBam(in_sam)
    sortBam(convertToBam.out) 
}