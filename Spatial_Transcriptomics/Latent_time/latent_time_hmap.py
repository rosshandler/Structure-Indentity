#############
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

seurat_prediction_colours = {
"Mitotic RG" : "#005dd8",
"Cycling RG" : "#4aadd6",
"Differentiating RG" : "#6083d1",
"Apical RG" : "#02d7ec",
"Glycolytic RG" : "#0086cb",
"High trancription RG" : "#4aac8d",
"Chp" : "#935de6",
"Cortical hem" : "#9e4d70",
"Cycling IPCs" : "#d74897",
"IPCs"   : "#ff8ba6",
"IN":"#228B22",
"Mature IN" : "#90b600",
"CR cells" : "#00ff00",
"Migrating excitatory neurons" : "#a6023e",
"UL neurons"   : "#fa2274",
"DL neurons"   : "#9e5d56",
"Mature excitatory neurons":"#d3b000",
}

genes_list = {
"PAX6",
"PTN",
"SOX2",
"HMGA2",
"SNTG1",
"FAM107A",
"LIFR",
"MOXD1",
"HOPX",
"TNC",
"EOMES",
"DLX6-AS1",
"CUX2",
"HTR2C",
"TTR",
"SLC17A7",
"SOX5",
"TBR1",
"PDZRN3",
"DLG2",
"GRIA2",
"NECAB1",
"TLE4",
"CACNA2D3",
"NR4A2",
"GAD2",
"GRIN2A",
"GRM7",
"ERBB4",
"RELN",
"GAD1",
"ETV1",
"NRF1",
"RORB",
"PTPRZ1",
"SATB2",
"BCL11B",
"NKX2-1",
"GRM1",
"CAMK2A",
"FOXP2",
"GRIK1",
} 

os.chdir('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')

###############################################

adata1 = sc.read('normcounts_postQC_slide1.tab')

file = open('cellnames_postQC_slide1.txt', 'r')
cells  = file.read().splitlines()

file = open('genenames_postQC_slide1.txt', 'r')
genes   = file.read().splitlines()

adata1.obs = pd.read_csv('metadata_postQC_slide1.csv', sep =',', low_memory=False)
adata1.obs_names = cells
adata1.var_names = genes

adata1.uns['seurat_prediction_colors'] = [seurat_prediction_colours[i] for i in sorted(np.unique(adata1.obs['seurat_prediction']))]

sc.pp.highly_variable_genes(adata1, n_top_genes=52)

hvgs1 = adata1.var_names[adata1.var['highly_variable']]

hvgs1 = list(filter(lambda x:'RPL' not in x, hvgs1))

scv.pl.heatmap(adata1, var_names=genes_list, sortby='latent_time', col_color='seurat_prediction', n_convolve=200, yticklabels=True)
hmap_ctrl = scv.pl.heatmap(adata1, var_names=genes_list, sortby='latent_time', col_color='seurat_prediction', n_convolve=200, yticklabels=True, show=False)
genes_ordered = hmap_ctrl.data.index.values

adata6 = sc.read('normcounts_postQC_slide6.tab')

file = open('cellnames_postQC_slide6.txt', 'r')
cells  = file.read().splitlines()

file = open('genenames_postQC_slide6.txt', 'r')
genes   = file.read().splitlines()

adata6.obs = pd.read_csv('metadata_postQC_slide6.csv', sep =',', low_memory=False)
adata6.obs_names = cells
adata6.var_names = genes

adata6.uns['seurat_prediction_colors'] = [seurat_prediction_colours[i] for i in sorted(np.unique(adata6.obs['seurat_prediction']))]

sc.pp.highly_variable_genes(adata6, n_top_genes=50)

hvgs6 = adata6.var_names[adata6.var['highly_variable']]

hvgs6 = list(filter(lambda x:'RPL' not in x, hvgs6))

scv.pl.heatmap(adata6, var_names=genes_list, sortby='latent_time', col_color='seurat_prediction', n_convolve=200, yticklabels=True)

scv.pl.heatmap(adata6, var_names=genes_ordered, sortby='latent_time', col_color='seurat_prediction', n_convolve=200, yticklabels=True, sort=False)

## Both together
adata=adata1.concatenate(adata6)

hvgs=list(set(hvgs1).union(hvgs6))

scv.pl.heatmap(adata, var_names=genes_list, sortby='latent_time', col_color='batch', n_convolve=300, yticklabels=True)
###

adata6 = sc.read('normcounts_postQC_slide6_A2.tab')

file = open('cellnames_postQC_slide6_A2.txt', 'r')
cells  = file.read().splitlines()

file = open('genenames_postQC_slide6_A2.txt', 'r')
genes   = file.read().splitlines()

adata6.obs = pd.read_csv('metadata_postQC_slide6_A2.csv', sep =',', low_memory=False)
adata6.obs_names = cells
adata6.var_names = genes

adata6.uns['seurat_prediction_colors'] = [seurat_prediction_colours[i] for i in sorted(np.unique(adata6.obs['seurat_prediction']))]

sc.pp.highly_variable_genes(adata6, n_top_genes=50)

hvgs6 = adata6.var_names[adata6.var['highly_variable']]

hvgs6 = list(filter(lambda x:'RPL' not in x, hvgs6))

scv.pl.heatmap(adata6, var_names=genes_list, sortby='latent_time', col_color='seurat_prediction', n_convolve=100, yticklabels=True)

