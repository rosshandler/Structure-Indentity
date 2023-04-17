import os
import numpy as np
import pandas as pd
import scanpy as sc
import sctour as sct
import matplotlib.pyplot as plt
os.chdir('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour')

adata = sc.read('raw_counts.mtx')

file = open('cells.txt', 'r')
cells  = file.read().splitlines()

file = open('genes.txt', 'r')
genes   = file.read().splitlines()

adata.obs =  pd.read_csv('metadata.tab',sep="\t")
adata.obs_names = cells
adata.var_names = genes

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
adata.obs

sc.pp.filter_genes(adata, min_cells=20)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)

tnode = sct.train.Trainer(adata, loss_mode='nb')
tnode.train()

adata.obs['ptime'] = tnode.get_time()

mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.2, alpha_predz=0.8)
adata.obsm['X_TNODE'] = mix_zs

adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15)
sc.tl.umap(adata, min_dist=0.1)
sc.tl.leiden(adata)
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])

fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))
sc.pl.umap(adata, color='leiden', size=20, ax=axs[0], legend_loc='on data', show=False)
sc.pl.umap(adata, color='ptime', size=20, ax=axs[1], show=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='leiden', ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2)
plt.show()