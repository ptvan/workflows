#!/usr/bin/bash

# map reads to reference genome
bwa mem -SP5M -t8 hg38.fast input_R1.fastq input_R2.fastq | \
   samtools view -bhS - > output.bam

# filter and sort mapped reads
samtools view -h output.bam | \
    pairtools parse -c hg38.mainonly.chrom.size -o parsed.pairsam.gz

pairtools sort --nproc 8 -o sorted.pairsam.gz parsed.pairsam.gz

# process pairsams: dedup and select
pairtools dedup --mark-dup -o deduped.pairsam.gz sorted.pairsam.gz
pairtools select \
    '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
    -o filtered.pairsam.gz deduped.pairsam.gz

# split pairsam and create pairs file
pairtools split --output-pairs output.pairs.gz filtered.pairsam.gz

# index pairs file
# NOTE: this operation is in-place and does not produce a new file
pairix -f output.pairs.gz

# binning into a .cool file
cooler cload pairix hg38.mainonly.chrom.size:500000 output.pairs.gz output.cool

# normalize .cool
# NOTE: in-place operation
cooler balance output.cool

# aggregate .cool into .mcool
# NOTE: in-place operation
cooler zoomify output.cool

# at this point the .mcool file can be view in HiGlass