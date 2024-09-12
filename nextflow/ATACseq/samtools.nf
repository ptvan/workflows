nextflow.enable.dsl=2

process SORTBAM {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(unsorted_bam_ch)

    output:
      tuple val(sample), path("${unsorted_bam_ch.baseName}.bam"), emit: sorted_bam_ch

    script:
    """
    samtools sort -O bam ${unsorted_bam_ch} -o ${unsorted_bam_ch.baseName}.tmp.bam
    samtools index ${unsorted_bam_ch.baseName}.tmp.bam -o ${unsorted_bam_ch.baseName}.bai
    mv ${unsorted_bam_ch.baseName}.tmp.bam ${unsorted_bam_ch.baseName}
    """    
}

process REMOVEMITOREADS{
publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(bam_mito_ch)

    output:
      tuple val(sample), path("${bam_mito_ch.baseName}.rmChrM.bam"), emit:bam_rmChrM_ch

    script:
    """
    samtools view -h ${bam_mito_ch} | grep -v chrM | samtools sort -O bam -o ${bam_mito_ch.baseName}.rmChrM.bam -T .
    samtools index ${bam_mito_ch.baseName}.rmChrM.bam -o ${bam_mito_ch.baseName}.rmChrM.bai
    """
}

process ADDREADGROUPS{
publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    
    // TO-DO : CUSTOMIZE EACH @RG FOR EACH SAMPLE !!!

    input:
      tuple val(sample), path(bam_noRG_ch)

    output:
      tuple val(sample), path("${bam_noRG_ch.baseName}.RGadded.bam"), emit: bam_RGadded_ch

    script:
    """
    samtools addreplacerg -r "@RG\tID:RG1\tSM:${sample}\tPL:Illumina\tLB:Library.fa" -o ${bam_noRG_ch.baseName}.RGadded.bam ${bam_noRG_ch.baseName}.bam
    samtools index ${bam_noRG_ch.baseName}.RGadded.bam -o ${bam_noRG_ch.baseName}.RGadded.bai
    """
}