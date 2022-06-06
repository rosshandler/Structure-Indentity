import os
import bbknn
import random
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

os.chdir('/data1/ivanir/Ilaria2021/data')

adata = sc.read('normalised_counts_qc_doublet_removal.tab')

file = open('barcodes_qc_doublet_removal.txt', 'r')
cells  = file.read().splitlines()

file = open('gene_ids_loom_filtered.txt', 'r')
genes   = file.read().splitlines()

adata.var_names_make_unique()

adata.obs = pd.read_csv('metadata_qc_doublet_removal.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

sc.pp.highly_variable_genes(adata, n_top_genes=250)

sc.pp.pca(adata, n_comps=30)

pca = pd.read_csv('pca_corrected_seurat.csv', sep =',', low_memory=False)
pca = pd.DataFrame(pca)

adata.obsm['X_pca'] = pca.values
adata.obsm['X_pca'] = adata.obsm['X_pca'].astype(float)

#umap = pd.read_csv('umap_seurat.csv', sep =',', low_memory=False)
#umap = pd.DataFrame(umap)
#adata.obsm['X_umap'] = umap.values
#adata.obsm['X_umap'] = adata.obsm['X_umap'].astype(float)

bbknn.bbknn(adata, batch_key = 'sequencing.round')

sc.tl.leiden(adata, resolution=1)

sc.tl.umap(adata, min_dist=.5)
np.savetxt('umap_layout_batch_corrected.csv', adata.obsm['X_umap'], delimiter=',')

cell_cycle_genes = [x.strip() for x in open('regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.obs.to_csv('structure-identity_clustering_batch_corrected.csv')

adata.obs[['scDblFinder.originAmbiguous']] = adata.obs[['scDblFinder.originAmbiguous']].astype('string')

sc.write('structure-identity_batch_corrected', adata)
