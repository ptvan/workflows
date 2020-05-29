# create coverage (in the form of a bigWig/bedGraph) file from a BAM file
bamCoverage -b reads.bam -o coverage.bw --numberOfProcessors 8

# scale mouse and fly GC contents and plot heatmaps for comparison
computeMatrix scale-regions -S mouse_GCcontent.bw -R RefSeq_genes_uniq.bed  -m 10000 -b 3000 -a 3000 -out mouse_matrix.tab.gz
plotHeatmap -m mouse_matrix.tab.gz -out scaledmouseGC.png --colorMap YlGnBu --regionsLabel 'mouse genes'

computeMatrix scale-regions -S fly_GCcontent.bw -R RefSeq_genes_uniq.bed -m 3000 -b 1000 -a 1000 -out fly_matrix.tab.gz
plotHeatmap -m fly_matrix.tab.gz -out scaledflyGC.png --colorMap YlGnBu --regionsLabel 'fly genes'