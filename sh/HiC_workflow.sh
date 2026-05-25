#!/usr/bin/bash

# map reads to reference genome
bwa mem -SP5M -t8 hg38.fast input_R1.fastq input_R2.fastq | \
   samtools view -bhS - > output.bam

# filter and sort mapped reads
pairtools parse -o output.parsed.pairsam.gz -c hg38.fa.sizes \
  --drop-sam --drop-seq --output-stats output.stats \
  --assembly hg38 --no-flip \
  --add-columns mapq \
  --walks-policy mask \
  output.bam

pairtools sort --nproc 8 -o output.sorted.pairs.gz output.parsed.pairs.gz

# deduplicate
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

# select
pairtools select "mapq1>0 and mapq2>0" output.nodups.pairs.gz -o test.nodups.UU.pairs.gz

# generate stats
pairtools stats output.sorted.pairs.gz -o output.stats

# create report from stats file
multiqc output.stats

# binning .pairs file into a .cool file
cooler cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    --assembly hg38 \
    ~/.local/share/genomes/hg38/hg38.fa.sizes:1000000 \ 
    output.sorted.pairs.gz \
    output.hg38.1000000.cool

# normalize .cool, aggregate into .mcool, perform balancing on each zoom level
# NOTE: this is an in-place operation and does not generate a new file
cooler zoomify \
    --nproc 5 \
    --out output.hg38.1000000.mcool \
    --resolutions 1000000,2000000 \
    --balance \
    output.hg38.1000000.cool

# at this point the .mcool file can be viewed in HiGlass or used to generate static figures
