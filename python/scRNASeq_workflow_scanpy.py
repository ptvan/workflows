import numpy as np
import pandas as pd
import scanpy as sc

# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write

# set params
sc.settings.verbosity = 3             
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = 'write/pbmc3k.h5ad'

# load data
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata.var_names_make_unique()

# plot genes that show highest expression
sc.pl.highest_expr_genes(adata, n_top=20, )

# flag cells and genes to filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# calculate metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# violin plots
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# scatter plots
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# perform filtering
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# normalize & take log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# identify highly-variable genes (HVG)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# plot HVGs
sc.pl.highly_variable_genes(adata)

# run PCA & plot results
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')

# inspect variance ratio
sc.pl.pca_variance_ratio(adata, log=True)

# save results
adata.write(results_file)

# compute neighbors graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# run UMAP
sc.tl.umap(adata)

# if there are disconnected clusters and other connectivity violations
#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)  
#sc.tl.umap(adata, init_pos='paga')

# plot UMAP results, w/ and w/o correction
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)

# cluster neighborhood graph
sc.tl.leiden(adata)

# plot clusters
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])

# find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
