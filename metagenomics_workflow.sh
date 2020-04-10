#!/usr/bin/bash

# get some COVID19 sequences from ViPR
# https://www.viprbrc.org/brc/home.spg?decorator=corona_ncov

# perform a multiple sequence alignment

# using MAFFT, outputting a .phy file
mafft --thread 8 --reorder --treeout COVID19multi.fa > COVID19multi.phy

# using ClustalW, outputting a .aln file
clustalw COVID19multi.fa

# construct evolutionary tree using IQ-TREE
# IQ-TREE takes either a .cf or .phy file 
# by default IQ-TREE will try to find the correct substitution model
# iqtree -nt 8 -s COVID19multi.phy

# looking at multiDNA.phy.iqtree we see that the best model is SYM+R5
# so on subsequent runs we can specify the model explicitly to save time
iqtree -nt 8 -s COVID19multi.phy -m HKY+F

# IQ-TREE output a .iqtree file (report) and .treefile file (NEWICK tree)



