import os
import numpy as np
import pandas as pd
import scanpy as sc

os.chdir('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy')

adata = sc.read('normalised_counts_qc.tab')
adata = adata.transpose()

file = open('cells_qc.txt', 'r')
cells  = file.read().splitlines()

file = open('genes_qc.txt', 'r')
genes   = file.read().splitlines()

adata.obs = pd.read_csv('cell_metadata_qc.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

adata.var_names_make_unique()

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.highly_variable_genes(adata, n_top_genes=600)

sc.pp.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(adata, resolution=1)

sc.tl.umap(adata, min_dist=.5)

sc.pl.umap(adata, color='phase')
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='condition')
sc.pl.umap(adata, color='seurat_prediction')

adata.write('structure-identity_pb')

#########
adata = sc.read('normalised_counts_ctrl.tab')
adata = adata.transpose()

file = open('cells_ctrl.txt', 'r')
cells  = file.read().splitlines()

file = open('genes_ctrl.txt', 'r')
genes   = file.read().splitlines()

adata.obs = pd.read_csv('cell_metadata_ctrl.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

adata.var_names_make_unique()

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.highly_variable_genes(adata, n_top_genes=600)

sc.pp.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(adata, resolution=1)

sc.tl.umap(adata, min_dist=.5)

sc.pl.umap(adata, color='phase')
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='condition')
sc.pl.umap(adata, color='seurat_prediction')

adata.write('structure-identity_ctrl')

#########
adata = sc.read('normalised_counts_diss.tab')
adata = adata.transpose()

file = open('cells_diss.txt', 'r')
cells  = file.read().splitlines()

file = open('genes_diss.txt', 'r')
genes   = file.read().splitlines()

adata.obs = pd.read_csv('cell_metadata_diss.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

adata.var_names_make_unique()

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.highly_variable_genes(adata, n_top_genes=600)

sc.pp.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(adata, resolution=1)

sc.tl.umap(adata, min_dist=.5)

sc.pl.umap(adata, color='phase')
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='condition')
sc.pl.umap(adata, color='seurat_prediction')

adata.write('structure-identity_diss')

#########
adata = sc.read('normalised_counts_emb.tab')
adata = adata.transpose()

file = open('cells_emb.txt', 'r')
cells  = file.read().splitlines()

file = open('genes_emb.txt', 'r')
genes   = file.read().splitlines()

adata.obs = pd.read_csv('cell_metadata_emb.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

adata.var_names_make_unique()

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.highly_variable_genes(adata, n_top_genes=600)

sc.pp.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(adata, resolution=1)

sc.tl.umap(adata, min_dist=.5)

sc.pl.umap(adata, color='phase')
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='condition')
sc.pl.umap(adata, color='seurat_prediction')

adata.write('structure-identity_emb')


######
sc.tl.diffmap(adata)
sc.pl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_diffmap', n_pcs=30)




sc.write('structure-identity_pb', adata)

adata.obs.to_csv('structure-identity_pb.csv')

np.savetxt('umap_layout.csv', adata.obsm['X_umap'], delimiter=',')
