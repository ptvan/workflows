#!/usr/bin/bash

# get example FASTA from MAFFT site
wget https://mafft.cbrc.jp/alignment/software/ex1.txt -O multiDNA.fa

# align sequences using MAFFT, output a .phy file
mafft --thread 8 --reorder --treeout multiDNA.fa > multiDNA.phy

# construct evolutionary tree using IQ-TREE
# IQ-TREE takes either a .cf or .phy file 
# by default IQ-TREE will try to find the correct substitution model
iqtree -nt 8 -s multiDNA.phy

# looking at multiDNA.phy.iqtree we see that the best model is SYM+R5
# so on subsequent runs we can specify the model explicitly to save time
iqtree -nt 8 -s multiDNA.phy -m SYM+R5

# IQ-TREE output a .iqtree file (report) and .treefile file (NEWICK tree)



