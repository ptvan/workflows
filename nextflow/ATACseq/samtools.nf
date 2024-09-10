process SORTBAM {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { unsorted_bam_ch.baseName }
    input:
      path unsorted_bam_ch

    output:
      path '*.sorted.ba*', emit: sortedbam_ch

    script:
    """
    samtools sort ${unsorted_bam_ch} > ${unsorted_bam_ch.baseName}.sorted.bam 
    samtools index ${unsorted_bam_ch.baseName}.sorted.bam -o ${unsorted_bam_ch.baseName}.sorted.bai
    """    
}

process REMOVEMITOREADS{
publishDir "${params.output}", mode:"copy", overwrite: true
    tag { bam_mito_ch.baseName }
    input:
      path bam_mito_ch

    output:
      path '*.rmChrM.ba*', emit: bam_rmChrM_ch

    script:
    """
    samtools view -h ${bam_mito_ch} | grep -v chrM | samtools sort -O bam -o ${bam_mito_ch.baseName}.rmChrM.bam -T .
    samtools index ${bam_mito_ch.baseName}.rmChrM.bam -o ${bam_mito_ch.baseName}.rmChrM.bai
    """

}