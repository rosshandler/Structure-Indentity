import os
import numpy as np
import pandas as pd
import scanpy as sc

os.chdir('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/Annotation')

adata = sc.read('normalised_counts.tab')
adata = adata.transpose()

file = open('cells.txt', 'r')
cells  = file.read().splitlines()

file = open('genes.txt', 'r')
genes   = file.read().splitlines()

adata.obs = pd.read_csv('metadata.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

sc.pp.highly_variable_genes(adata, n_top_genes=250)

sc.pp.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.leiden(adata, resolution=1)

sc.tl.umap(adata, min_dist=.5)

sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='day')
sc.pl.umap(adata, color='condition')
sc.pl.umap(adata, color='seurat_prediction')




sc.write('structure-identity_batch_corrected', adata)

adata.obs.to_csv('structure-identity_clustering_batch_corrected.csv')

np.savetxt('umap_layout_batch_corrected.csv', adata.obsm['X_umap'], delimiter=',')
