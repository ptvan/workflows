nextflow.enable.dsl=2

process FILTERBLACKLISTREGIONS {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(input_bam_ch)
      path(bed_file)
  
    output:      
      tuple val(sample), path("${input_bam_ch.baseName}.blacklistfiltered.bam"), emit: blacklistfiltered_ch

    script:
    """
    bedtools intersect -nonamecheck -v -abam ${input_bam_ch} -b ${bed_file} > ${input_bam_ch.baseName}.blacklistfiltered.bam
    """    
}

process BAMTOBED {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(input_bam_ch)
      
    output:
      path("${bed_ch.baseName}.bed"), emit: bed_ch

    script:
    """
    bedtools bamtobed -i ${input_bam_ch} > ${bed_ch.baseName}.bed
    """    
}

process BAMTOBEDPE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(input_bam_ch)
      
    output:
      path("${bed_ch.baseName}.bed"), emit: bed_ch


    script:
    """
    bedtools bamtobed -i ${input_bam_ch} -bedpe > ${bedpe_ch.baseName}.bedpe
    """    
}