nextflow.enable.dsl=2

process MARKBLACKLISTREGIONS {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { input_bam_ch.baseName }
    input:
      path input_bam_ch
      path bed_file

    output:
      path '*.blacklistfiltered.bam', emit: blacklisted_bam_ch

    script:
    """
    bedtools intersect -nonamecheck -v -abam ${input_bam_ch} -b ${bed_file} > ${input_bam_ch.baseName}.blacklisted.bam
    """    
}

process BAMTOBED {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { input_bam_ch.baseName }
    input:
      path input_bam_ch
      
    output:
      path '*.bed', emit: bed_ch

    script:
    """
    bedtools bamtobed -i ${input_bam_ch} > ${bed_ch.baseName}.bed
    """    
}

process BAMTOBEDPE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { input_bam_ch.baseName }
    input:
      path input_bam_ch
      
    output:
      path '*.bedpe', emit: bedpe_ch

    script:
    """
    bedtools bamtobed -i ${input_bam_ch} -bedpe > ${bedpe_ch.baseName}.bedpe
    """    
}