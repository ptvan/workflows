// based on https://datacarpentry.org/wrangling-genomics/
nextflow.enable.dsl=2

params.suffix = ".fastq.gz"
params.help = false

if (params.help || params.input == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run multiQC --input <> --output <>
    
    Required Arguments:
      --input        Folder containing all input data in FASTQ files (will traverse subdirectories)
      --output       Folder to place analysis outputs (named 'multiqc_report.html')
    Input Files:
      --suffix              Process all files ending with this string (default: .fastq.gz)
    """.stripIndent()
}

// Show help message if --help flag at runtime

process fastQC {
  input:
    tuple val(name), file(reads)

  output:
    file "*_fastqc.zip"

"""
#!/bin/bash
set -Eeuo pipefail
fastqc -t ${task.cpus} -o ./ ${reads}
ls -lahtr
"""

}

process multiQC {

  publishDir "${params.output}", mode: 'copy', overwrite: true
  
  input:
    file "*_fastqc.zip"

  output:
    file "multiqc_report.html"

"""
#!/bin/bash
set -Eeuo pipefail
multiqc .
ls -lahtr
"""

}

// Start the workflow
workflow {

    // Get the input files ending with BAM
    input_ch = Channel.fromPath(
        "${params.input}**${params.suffix}"
    ).map {
        it -> [it.name.replaceAll(/${params.suffix}/, ''), it]
    }

    // Run FastQC on the reads
    fastQC(
        input_ch
    )

    // Aggregate all results with MultiQC
    multiQC(
        fastQC.out.toSortedList()
    )

}