// based on https://datacarpentry.org/wrangling-genomics/
nextflow.enable.dsl=2

include { RUNFASTQC } from 'runFastQC'
include { TRIMADAPTERS } from 'trimAdapters'
include { INDEXREFERENCE } from 'indexReferenceGenome'


workflow wranglingGenomics {
    take:
        untrimmed_fastq
        ref
        output
    
    main:
        RUNFASTQC(
            untrimmed_fastq
        )
        TRIMADAPTERS(
            untrimmed_fastq
        )
        INDEXREFERENCE(
            ref
        )

}
