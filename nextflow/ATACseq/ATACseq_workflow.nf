nextflow.enable.dsl=2
params.suffix = "*_R{1,2}.trimmed.fastq.gz"
params.alignmentProgram = "/usr/local/bin/bowtie2"
params.alignmentParams = "--local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x "
params.referenceGenome = "$HOME/working/Databases/GCRh38_ATACseq"
params.genomeBlacklist = "$HOME/working/raw_data/hg38.blacklist.bed.gz"
params.nCPUs = 8
params.help = false

include { ALIGNTOREFERENCE } from './bowtie2.nf'
include { SORTBAM; REMOVEMITOREADS } from './samtools.nf'
include { REMOVEDUPLICATEREADS } from './picard.nf'
include { MARKBLACKLISTREGIONS } from './bedtools.nf'
include { RUNALIGNMENTSIEVE; RUNBAMCOVERAGE } from './deeptools.nf'
include { CALLPEAKS } from './MACS.nf'

if (params.help || params.input == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run align.nf --input <> --output <>
        
    NOTE: this pipeline requires GNU parallel 
    
    Required Arguments:
      --input        Path to .fasta.gz input files 
      --output       Path to store output files 
    
    """.stripIndent()
}

workflow {
    raw_reads = Channel.fromFilePairs("${params.input}/${params.suffix}")
    aligned_reads = ALIGNTOREFERENCE(raw_reads)
    sorted_reads = SORTBAM(aligned_reads)
    noMito_reads = REMOVEMITOREADS(sorted_reads)
}