#!/usr/bin/bash

# get example FASTA from MAFFT site
wget https://mafft.cbrc.jp/alignment/software/ex1.txt -O multiDNA.fa

# align sequences using MAFFT
mafft --thread 8 --reorder --treeout multiDNA.fa > multiDNA.phy


# construct evolutionary tree using IQ-TREE
# IQ-TREE takes either a .cf or .phy file 
# and output a .iqtree file (report) and .treefile file (NEWICK tree)
iqtree -nt 8 -redo -s example.phy
