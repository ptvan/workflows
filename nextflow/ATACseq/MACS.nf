nextflow.enable.dsl=2

process RANDSAMPLE {
  publishDir "${params.output}", mode:"copy", overwrite: true
  tag { sample }
  input:
    tuple val(sample), path(bam_ch)

  output:
    tuple val(sample), path("*.bed"), emit: randsampleoutput_ch
  
  script:
  """
  macs3 randsample \
  -i ${bam_ch} \
  -f BAMPE \
  -p 100 \
  -o ${sample}.bed

  """
}

process CALLPEAKS {
  publishDir "${params.output}", mode:"copy", overwrite: true
  tag { sample }
  input:
    tuple val(sample), path(bam_ch)

  output:
    tuple val(sample), path("MACS/*"), emit: callpeaksoutput_ch

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
  --cutoff-analysis -n ${sample} \
  -t ${bam_ch.baseName}.bed \
  --outdir ${params.output} 2> ${bam_ch.baseName}.macs3.log
  """
}
