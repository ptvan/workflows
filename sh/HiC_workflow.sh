#!/usr/bin/bash

# create the <GENOME>.sizes file required by ligation and binning steps
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.fa.sizes

# map reads to reference genome using BWA, disable pairing
bwa mem -SP5M -t8 hg38.fa input_R1.fastq input_R2.fastq | \
   samtools view -bhS - > output.bam

# detect ligation pairs from mapped reads
pairtools parse -o output.parsed.pairs.gz -c hg38.fa.sizes \
  --drop-sam --drop-seq --output-stats output.stats \
  --assembly hg38 --no-flip \
  --add-columns mapq \
  --walks-policy mask \
  output.bam

# sort pairs lexicographically
pairtools sort --nproc 8 -o output.sorted.pairs.gz output.parsed.pairs.gz

# deduplicate pairs
pairtools dedup \
    --max-mismatch 3 \
    --mark-dups \
    --output \
        >( pairtools split \
            --output-pairs output.nodups.pairs.gz \
            --output-sam output.nodups.bam \
         ) \
    --output-unmapped \
        >( pairtools split \
            --output-pairs output.unmapped.pairs.gz \
            --output-sam output.unmapped.bam \
         ) \
    --output-dups \
        >( pairtools split \
            --output-pairs output.dups.pairs.gz \
            --output-sam output.dups.bam \
         ) \
    --output-stats output.dedup.stats \
    output.sorted.pairs.gz

# filter pairs
# pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' output.nodups.pairs.gz -o output.nodups.UU.pairs.gz
pairtools select "mapq1>0 and mapq2>0" output.nodups.pairs.gz -o output.nodups.UU.pairs.gz

# generate stats
pairtools stats output.nodups.UU.pairs.gz -o output.stats

# create report from stats file
multiqc output.stats

# binning .pairs file into a .cool file
cooler cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    --assembly hg38 \
    hg38.fa.sizes:1000000 \ 
    output.nodups.UU.pairs.gz \
    output.hg38.1000000.cool

# aggregate .cool into .mcool, perform normalization on each zoom level
# NOTE: this is an in-place operation and does not generate a new file
cooler zoomify \
    --nproc 5 \
    --out output.hg38.1000000.mcool \
    --resolutions 1000000,2000000 \
    --balance \
    output.hg38.1000000.cool

# at this point the .mcool file can be viewed in HiGlass or used to generate static figures
