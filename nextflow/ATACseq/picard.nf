nextflow.enable.dsl=2
params.picardProgram = "~/working/packages/picard.jar"

process REMOVEDUPLICATEREADS {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { bam_dupes_ch.baseName }
    input:
      path bam_dupes_ch
   
    output:
      path '*.bam', emit: bam_nodupes_ch

    script:
    
    """
    ${params.picardProgram} MarkDuplicates \
    QUIET=true \ 
    INPUT=${bam_dupes_ch} \
    OUTPUT=${bam_dupes_ch}.dupesremoved.bam \
    METRICS_FILE=${bam_dupes_ch}.dup.metrics \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=.
    
    """
}