nextflow.enable.dsl=2
// TO-DO : HAVE LOGGING FOR EACH SAMPLE !!!
process CALLPEAKS {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { bam_ch.baseName }
    input:
      path bam_ch

    output:
      path 'MACS/*', emit: MACS_output_ch

    script:
    """
    macs3 callpeak \
    -f BEDPE \
    --nomodel \
    --shift -37 \
    --extsize 73 \
    -g 2862010578 \
    -B --broad \
    --keep-dup all \
    --cutoff-analysis -n ${bam_ch} \
    -t ${bam_ch.baseName}.bed \
    --outdir ${params.output} 2> ${bam_ch.baseName}.macs3.log
    """
}
