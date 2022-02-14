import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

Cluster_colours = {
  "Mic"   : "#935de6",
  "InMGE" : "#90b600",
  "PgG2M" : "#005dd8",
  "ExDp1" : "#d3b000",
  "PgS"   : "#5484ff",
  "InCGE" : "#00c571",
  "ExM-U" : "#fa2274",
  "oRG"   : "#02d7ec",
  "ExM"   : "#ff6500",
  "OPC"   : "#02acf4",
  "End"   : "#bd6200",
  "vRG"   : "#ff7bd8",
  "ExDp2" : "#9e5d56",
  "Per"   : "#ada7de",
  "ExN"   : "#a6023e",
  "IP"    : "#ff8ba6",
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

os.chdir('/data1/ivanir/HumaFetalNeocortex2019/data')

###############################################

adata1 = sc.read('normalised_counts_invivo.tab')

file = open('barcodes_invivo.txt', 'r')
cells  = file.read().splitlines()

file = open('gene_ids_invivo.txt', 'r')
genes   = file.read().splitlines()

adata1.obs = pd.read_csv('meta_ltime_invivo.csv', sep =',', low_memory=False)
adata1.obs_names = cells
adata1.var_names = genes

adata1.uns['Cluster_colors'] = [Cluster_colours[i] for i in sorted(np.unique(adata1.obs['Cluster']))]

sc.pp.highly_variable_genes(adata1, n_top_genes=60)

hvgs1 = adata1.var_names[adata1.var['highly_variable']]

hvgs1 = list(filter(lambda x:'RPL' not in x, hvgs1))

scv.pl.heatmap(adata1, var_names=genes_list, sortby='projected_latent_time', col_color='Cluster', n_convolve=500, yticklabels=True)

scv.pl.heatmap(adata1, var_names=hvgs1, sortby='projected_latent_time', col_color='Cluster', n_convolve=500, yticklabels=True)

def Union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list
  
genes_list_extended = Union(genes_list, hvgs1)
scv.pl.heatmap(adata1, var_names=genes_list_extended, sortby='projected_latent_time', col_color='Cluster', n_convolve=500, yticklabels=True)


###############################################

adata1 = sc.read('normalised_counts_invivo_to_morph_cluster1.tab')

file = open('barcodes_invivo_to_morph_cluster1.txt', 'r')
cells  = file.read().splitlines()

file = open('gene_ids_invivo_to_morph_cluster1.txt', 'r')
genes   = file.read().splitlines()

adata1.obs = pd.read_csv('meta_ltime_invivo_to_morph_cluster1.csv', sep =',', low_memory=False)
adata1.obs_names = cells
adata1.var_names = genes

adata1.uns['Cluster_colors'] = [Cluster_colours[i] for i in sorted(np.unique(adata1.obs['Cluster']))]

sc.pp.highly_variable_genes(adata1, n_top_genes=60)

hvgs1 = adata1.var_names[adata1.var['highly_variable']]

hvgs1 = list(filter(lambda x:'RPL' not in x, hvgs1))

scv.pl.heatmap(adata1, var_names=genes_list, sortby='projected_latent_time', col_color='Cluster', n_convolve=500, yticklabels=True)


###############################################

adata1 = sc.read('normalised_counts_invivo_to_morph_cluster2.tab')

file = open('barcodes_invivo_to_morph_cluster2.txt', 'r')
cells  = file.read().splitlines()

file = open('gene_ids_invivo_to_morph_cluster2.txt', 'r')
genes   = file.read().splitlines()

adata1.obs = pd.read_csv('meta_ltime_invivo_to_morph_cluster2.csv', sep =',', low_memory=False)
adata1.obs_names = cells
adata1.var_names = genes

adata1.uns['Cluster_colors'] = [Cluster_colours[i] for i in sorted(np.unique(adata1.obs['Cluster']))]

sc.pp.highly_variable_genes(adata1, n_top_genes=60)

hvgs1 = adata1.var_names[adata1.var['highly_variable']]

hvgs1 = list(filter(lambda x:'RPL' not in x, hvgs1))

scv.pl.heatmap(adata1, var_names=genes_list, sortby='projected_latent_time', col_color='Cluster', n_convolve=500, yticklabels=True)



