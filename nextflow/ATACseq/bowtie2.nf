nextflow.enable.dsl=2

process ALIGNTOREFERENCE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), file(raw_reads)
   
    output:
      path '*.bam', emit: unsorted_bam_ch

    script:
    def (read1, read2) = raw_reads
    """
    ${params.alignmentProgram} ${params.alignmentParams} ${params.referenceGenome} \
        -1 ${read1} -2 ${read2} -S ${sample}.bam
    
    """
}