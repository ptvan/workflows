nextflow.enable.dsl=2

params.suffix = "*_R{1,2}.trimmed.fastq.gz"
params.alignmentProgram = "/usr/local/bin/bowtie2"
params.alignmentParams = "--local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x "
params.referenceGenome = "$HOME/working/Databases/GCRh38_ATACseq"
params.help = false

if (params.help || params.input == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run align.nf --input <> --output <>
    
    Required Arguments:
      --input        Full path to .fasta.gz file to align
      --output       Folder to place BAM files
    
    NOTE: this pipeline requires GNU parallel 
    """.stripIndent()
}

process alignToReference {
    publishDir "${params.output}", mode:"copy", overwrite: true

    input:
      tuple val(sample), file(raw_reads)
   
    output:
      path '*.bam', emit: unsorted_bam_ch

    script:
    def (read1, read2) = raw_reads
    """
    ${params.alignmentProgram} ${params.alignmentParams} ${params.referenceGenome} \
        -1 ${read1} -2 ${read2} | samtools view -bS - > ${sample}.bam
    
    """
}

process sortBam {
    publishDir "${params.output}", mode:"copy", overwrite: true

    input:
      path unsorted_bam_ch

    output:
      path '*.sorted.bam', emit: sortedbam_ch

    script:
    """
    ls *.bam | parallel -j 8 'samtools sort {} > {.}.sorted.bam'
    """    
}

workflow {
    raw_reads = Channel.fromFilePairs("${params.input}/${params.suffix}")
    aligned_reads = alignToReference(raw_reads)
    sorted_reads = sortBam(aligned_reads) 
}