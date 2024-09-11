nextflow.enable.dsl=2

process RUNALIGNMENTSIEVE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { bam_ch.baseName }
    input:
      path bam_ch

    output:
      path '*.shifted.ba*', emit: shiftedbam_ch

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
    tag { bam_ch.baseName }
    input:
      path bam_ch

    output:
      path '*.bw', emit: shiftedbam_ch

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