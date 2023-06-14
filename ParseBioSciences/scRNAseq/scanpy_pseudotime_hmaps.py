import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

os.chdir('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy')

from colour import Color
blue = Color("blue")
colors = list(blue.range_to(Color("green"),10))

# colors is now a list of length 10
# Containing: 
#[<Color blue>, <Color #0036f1>, <Color #0065e3>, <Color #008ed5>, <Color #00b0c7>, <Color #00b8a4>, <Color #00aa72>, <Color #009c45>, <Color #008e20>, <Color green>]

leiden_annot_colours = {
"Differentiating RG" : "#d54d92",
"Mixed Identity RG/Neu 1" : "#5cc556",
"Mixed Identity RG/Neu 2" : "#b254bf",
"Chp 1" : "#9cb735",
"Glicolytic Neuronal" : "#6c69ca",
"RG 1" : "#cda937",
"IPCs" : "#5d8fcb",
"Migrating Excitatory Neurons" : "#dc5b31",
"Mitotic RG 1" : "#42bdc0",
"Mixed Identity RG/Neu 3" : "#dd4663",
"UL Neurons" : "#56993f",
"DL Neurons 1" : "#c78acc",
"Mixed Identity RG/Neu 4" : "#61bf8c",
"RG 2" : "#b63f37",
"Chp 2" : "#407a48",
"Inhibitory Neurons" : "#9e4a6b",
"DL Neurons 2" : "#b5ae68",
"CR Cells" : "#e18880",
"Mitotic RG 2" : "#757327",
"Mixed Identity Chp/Neu" : "#d48a3c",
"Chemochine Signaling" : "#9c5e32",
}

day_colours = {
"45 day" : "#0036f1",
"55 day" : "#00b8a4",
"70 day" : "#008e20",
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
#"HOPX",
"TNC",
"EOMES",
#"DLX6-AS1",
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
#"PTPRZ1",
"SATB2",
"BCL11B",
#"NKX2-1",
"GRM1",
"CAMK2A",
"FOXP2",
"GRIK1",
} 

adata_ctrl = sc.read('normalised_counts_ctrl.tab').T

file = open('cells_ctrl.txt', 'r')
cells  = file.read().splitlines()

file = open('genes_ctrl.txt', 'r')
genes   = file.read().splitlines()

adata_ctrl.obs = pd.read_csv('cell_metadata_ctrl.tab', sep ='\t', low_memory=False)
adata_ctrl.obs_names = cells
adata_ctrl.var_names = genes
adata_ctrl.var_names_make_unique() 

adata_ctrl.obs['day'].value_counts()
#55 day    2815
#45 day    1789
#70 day     900
#Name: day, dtype: int64

adata_ctrl.uns['day'] = [day_colours[i] for i in sorted(np.unique(adata_ctrl.obs['day']))]
adata_ctrl.uns['leiden_annot'] = [leiden_annot_colours[i] for i in sorted(np.unique(adata_ctrl.obs['leiden_annot']))]

sc.tl.pca(adata_ctrl, svd_solver='arpack', random_state=1)
sc.pp.neighbors(adata_ctrl, n_neighbors=30)
sc.tl.umap(adata_ctrl)

sc.pl.umap(adata_ctrl, color='day', palette=day_colours)
sc.pl.umap(adata_ctrl, color='leiden_annot',palette=leiden_annot_colours)

scv.pl.heatmap(adata_ctrl, var_names=genes_list, sortby='pt_monocle3', col_color='day', palette=day_colours, n_convolve=200, yticklabels=True)
scv.pl.heatmap(adata_ctrl, var_names=genes_list, sortby='pt_monocle3', col_color='leiden_annot', palette=leiden_annot_colours, n_convolve=200, yticklabels=True)

hmap_ctrl = scv.pl.heatmap(adata_ctrl, var_names=genes_list, sortby='pt_monocle3', col_color='day', n_convolve=200, yticklabels=True, show=False)
genes_ordered = hmap_ctrl.data.index.values

###
day_colours = {
"48 day" : "#0036f1",
"55 day" : "#00b8a4",
"70 day" : "#008e20",
}

adata_diss = sc.read('normalised_counts_diss.tab').T

file = open('cells_diss.txt', 'r')
cells  = file.read().splitlines()

file = open('genes_diss.txt', 'r')
genes   = file.read().splitlines()

adata_diss.obs = pd.read_csv('cell_metadata_diss.tab', sep ='\t', low_memory=False)
adata_diss.obs_names = cells
adata_diss.var_names = genes
adata_diss.var_names_make_unique() 

adata_diss.obs['day'].value_counts()

sc.tl.pca(adata_diss, svd_solver='arpack', random_state=1)
sc.pp.neighbors(adata_diss, n_neighbors=30)
sc.tl.umap(adata_diss)
adata_diss.uns['day'] = [day_colours[i] for i in sorted(np.unique(adata_diss['day']))]

sc.pl.umap(adata_diss, color='day', palette=day_colours)

scv.pl.heatmap(adata_diss, var_names=genes_ordered, sortby='pt_monocle3', col_color='day', n_convolve=200, palette=day_colours,yticklabels=True, sort=False)
scv.pl.heatmap(adata_diss, var_names=genes_list, sortby='pt_monocle3', col_color='day', palette=day_colours,n_convolve=200, yticklabels=True)









import pandas as pd
import scanpy as sc
import numpy as np
import bbknn  
import os

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
sc.pp.highly_variable_genes(adata, n_top_genes=650)

sc.tl.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=1)
adata.obs = adata.obs.rename({'leiden': 'leiden_pca'}, axis='columns')

bbknn.ridge_regression(adata, batch_key=['batch'],confounder_key=['leiden_pca'])
sc.tl.pca(adata, n_comps=30)
bbknn.bbknn(adata)

#sc.tl.umap(adata)
#sc.pl.umap(adata, color='seurat_prediction')

sc.tl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap', n_pcs=30)
sc.tl.leiden(adata, resolution=1)
adata.obs = adata.obs.rename({'leiden': 'leiden_dmap'}, axis='columns')

adata.obs["seurat_prediction"] = adata.obs["seurat_prediction"].astype("category")
sc.tl.paga(adata, groups='seurat_prediction')
sc.pl.paga(adata)

sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color='seurat_prediction')
#sc.pl.umap(adata, color='leiden_dmap')
#sc.pl.umap(adata, color='leiden_pca')
#sc.pl.umap(adata, color='condition')
#sc.pl.umap(adata, color='phase')

adata.write('structure-identity_pb.hda5')

###############
###############
import pandas as pd
import scanpy as sc
import numpy as np
import bbknn  
import os

os.chdir('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy')

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
sc.pp.highly_variable_genes(adata, n_top_genes=650)

sc.tl.pca(adata, n_comps=30)

cell_cycle_genes = [x.strip() for x in open('/data1/ivanir/Ilaria2021/data/regev_lab_cell_cycle_genes.txt')]
s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=1)
adata.obs = adata.obs.rename({'leiden': 'leiden_pca'}, axis='columns')

bbknn.ridge_regression(adata, batch_key=['batch'],confounder_key=['leiden_pca'])
sc.tl.pca(adata, n_comps=30)
bbknn.bbknn(adata)

#sc.tl.umap(adata)
#sc.pl.umap(adata, color='seurat_prediction')

sc.tl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap', n_pcs=30)
sc.tl.leiden(adata, resolution=1)
adata.obs = adata.obs.rename({'leiden': 'leiden_dmap'}, axis='columns')

adata.obs["seurat_prediction"] = adata.obs["seurat_prediction"].astype("category")
sc.tl.paga(adata, groups='seurat_prediction')
sc.pl.paga(adata)

sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color='seurat_prediction')
#sc.pl.umap(adata, color='leiden_dmap')
#sc.pl.umap(adata, color='leiden_pca')
#sc.pl.umap(adata, color='condition')
#sc.pl.umap(adata, color='phase')

adata.write('structure-identity_pb_ctrl.hda5')

adata.obs.to_csv('structure-identity_pb_ctrl.csv')

np.savetxt('umap_layout_ctrl.csv', adata.obsm['X_umap'], delimiter=',')



