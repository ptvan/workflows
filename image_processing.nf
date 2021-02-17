#!/usr/bin/env nextflow
// convert the images in the CMU Face Images data set into JPEG
// original: https://archive.ics.uci.edu/ml/datasets/cmu+face+images

results_path = "~/working/aws/results"

params.input = "~/working/aws/faces/an2i/*pgm"
files = Channel.fromPath(params.input)

process convert {
    publishDir "$results_path/"
    input:
    file x from files

    output:
    set x, file("${x.baseName}.jpg") into converted

    script:
    """
    convert ${x} ${x.baseName}.jpg
    """
}


