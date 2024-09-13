nextflow.enable.dsl=2

process SORTBAM {
  publishDir "${params.output}", mode:"copy", overwrite: true
  tag { sample }
  input:
    tuple val(sample), path(unsorted_bam_ch)
    tuple val(sample), path(unsorted_bai_ch)

  output:
    tuple val(sample), path("${unsorted_bam_ch.baseName}.bam"), emit: sorted_bam_ch
    tuple val(sample), path("${unsorted_bai_ch.baseName}.bai"), emit: sorted_bai_ch

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
    tuple val(sample), path(bai_mito_ch)

  output:
    tuple val(sample), path("${bam_mito_ch.baseName}.noChrM.bam"), emit:bam_noChrM_ch
    tuple val(sample), path("${bai_mito_ch.baseName}.noChrM.bai"), emit:bai_noChrM_ch

  script:
  """
  samtools view -h ${bam_mito_ch} | grep -v chrM | samtools sort -O bam -o ${bam_mito_ch.baseName}.noChrM.bam -T .
  samtools index ${bam_mito_ch.baseName}.noChrM.bam -o ${bam_mito_ch.baseName}.noChrM.bai
  """
}

process ADDREADGROUPS{
  publishDir "${params.output}", mode:"copy", overwrite: true
  tag { sample }
  input:
    tuple val(sample), path(bam_noRG_ch)
    tuple val(sample), path(bai_noRG_ch)

  output:
    tuple val(sample), path("${bam_noRG_ch.baseName}.RGadded.bam"), emit: bam_RGadded_ch
    tuple val(sample), path("${bai_noRG_ch.baseName}.RGadded.bai"), emit: bai_RGadded_ch

  script:
  """
  samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}\tPL:Illumina\tLB:Library.fa" -o ${bam_noRG_ch.baseName}.RGadded.bam ${bam_noRG_ch.baseName}.bam
  samtools index ${bam_noRG_ch.baseName}.RGadded.bam -o ${bai_noRG_ch.baseName}.RGadded.bai
  """
}