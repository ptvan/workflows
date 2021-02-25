#!/usr/bin/env nextflow
 
params.in = "$HOME/working/datasets/parsnp-examples/mers_virus/genomes/"
params.ref = "$HOME/working/datasets/parsnp-examples/mers_virus/ref/England1.gbk"

input_files = file(params.in)
db = file(params.ref)

process runParsnp {
    
    input:
    path '*.fna' from params.in
 
    output:
    stdout results
 
    """
    parsnp -g $db -d $input_files -c
    """

}
