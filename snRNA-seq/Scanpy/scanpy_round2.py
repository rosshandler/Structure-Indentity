import os
import bbknn
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

os.chdir('/data1/ivanir/Ilaria2021/data')

adata = sc.read('structure-identity_batch_corrected')

file = open('cell_barcodes_scvelo.txt', 'r')
cells  = file.read().splitlines()

keep_cells = np.in1d(adata.obs_names,cells)

adata = adata[keep_cells]

sc.tl.umap(adata, min_dist=.5)

np.savetxt('umap_layout_batch_corrected_contamination_cleaned.csv', adata.obsm['X_umap'], delimiter=',')
