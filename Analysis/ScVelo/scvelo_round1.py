######
import os
import numpy as np
import loompy as lp
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad

import scipy.sparse as sparse
import scipy.io as sio
import scipy.stats as stats


os.chdir('/data1/ivanir/Ilaria2021/data')

lp.combine(
["/data2/mlancast/cruk/SLX-19865/fastq/B6Fast/velocyto/B6Fast.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B6Kit/velocyto/B6Kit.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B7diss/velocyto/B7diss.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B7Fast/velocyto/B7Fast.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B7Kit/velocyto/B7Kit.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B8diss/velocyto/B8diss.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B8Kit/velocyto/B8Kit.loom",
 "/data2/mlancast/cruk/SLX-19865/fastq/B8Unembed/velocyto/B8Unembed.loom",
 "/data2/mlancast/cruk/SLX-20646/B6diss/velocyto/B6diss.loom",
 "/data2/mlancast/cruk/SLX-20646/B6Unembed/velocyto/B6Unembed.loom",
 "/data2/mlancast/cruk/SLX-20646/B8Fast/velocyto/B8Fast.loom",
 "/data2/mlancast/cruk/SLX-20646/B8Unembed/velocyto/B8Unembed.loom"],
 "structure-identity.loom")
 
adata = sc.read('structure-identity.loom', sparse=True)

sio.mmwrite("unspliced_matrix.mtx",adata.layers['unspliced'].T, field='integer')
sio.mmwrite("spliced_matrix.mtx",adata.layers['spliced'].T, field='integer')

np.savetxt('barcodes_loom.csv', adata.obs_names.values,fmt='%5s', delimiter=',')
np.savetxt('gene_ids_loom.csv', adata.var_names.values,fmt='%5s', delimiter=',')

pd.DataFrame.sparse.from_spmatrix(adata.layers['unspliced'].T, index=adata.var_names, columns=adata.obs_names).to_csv('unspliced_raw_counts.csv')
pd.DataFrame.sparse.from_spmatrix(adata.layers['spliced'].T, index=adata.var_names, columns=adata.obs_names).to_csv('spliced_raw_counts.csv')

file = open('barcodes_qc_doublet_removal_loom.txt', 'r')
barcodes = file.read().splitlines()

adata.obs_names_make_unique()
adata.var_names_make_unique()

adata = adata[barcodes]
adata.obs = pd.read_csv('meta_annotated_updated.csv', low_memory=False)

adata.obs_names_make_unique()
adata.var_names_make_unique()

scv.pl.proportions(adata)

#scv.pp.pca(adata, n_comps=30)
#pca = pd.read_csv('pca_corrected_seurat.csv', sep =',', low_memory=False)
#pca = pd.DataFrame(pca)
#adata.obsm['X_pca'] = pca.values
#adata.obsm['X_pca'] = adata.obsm['X_pca'].astype(float)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=5000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.umap(adata)

umap = pd.read_csv('umap_layout_batch_corrected.csv', sep =',', low_memory=False)
umap = pd.DataFrame(umap)
adata.obsm['X_umap'] = umap.values
adata.obsm['X_umap'] = adata.obsm['X_umap'].astype(float)

adata.obs[['scDblFinder.originAmbiguous']] = adata.obs[['scDblFinder.originAmbiguous']].astype('string')

scv.tl.recover_dynamics(adata, n_jobs=60)
scv.tl.velocity(adata, mode='dynamical')#has to be done in scvelo-p3.6
adata.write('./write/scvelo_analysis_dynamical.h5ad')
