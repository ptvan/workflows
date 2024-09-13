nextflow.enable.dsl=2

process RUNALIGNMENTSIEVE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(bam_ch)

    output:
      tuple val(sample), path("${bam_noRG_ch.baseName}.shifted.bam"), emit: shiftedbam_ch

    script:
    """
    alignmentSieve \
    --verbose \
    --ATACshift \
    --blackListFileName ${params.genomeBlacklist} \
    --bam ${bam_ch}.bam \
    -o ${bam_ch}.shifted.bam
    
    samtools index ${bam_ch.baseName}.shifted.bam -o ${bam_ch.baseName}.shifted.bai
    """
}

process RUNBAMCOVERAGE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(bam_ch)

    output:
      tuple val(sample), path("*.bw"), emit: coveragebigwig_ch

    script:
    """
    bamCoverage \
    --numberOfProcessors ${params.nCPUs} \
    --binSize 10 \
    --normalizeUsing BPM \
    --effectiveGenomeSize 2862010578 \
    --bam ${bam_ch}.bam \
    -o ${bam_ch}_coverage_BPM.bw
        
    """
}