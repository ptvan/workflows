nextflow.enable.dsl=2
params.picardMarkDuplicates = "java -jar ~/working/packages/picard.jar MarkDuplicates"

process REMOVEDUPLICATEREADS {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
     tuple val(sample), path(bam_dupes_ch)
   
    output:
      tuple val(sample), path("${bam_dupes_ch.baseName}.rmDuplicates.bam"), emit: bam_nodupes_ch

    script:
    
    """
    ${params.picardMarkDuplicates} QUIET=true INPUT=${bam_dupes_ch} OUTPUT=${bam_dupes_ch.baseName}.rmDuplicates.bam METRICS_FILE=${bam_dupes_ch}.duplicates.metrics REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
    
    """
}