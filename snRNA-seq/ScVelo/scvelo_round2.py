######
import os
import numpy as np
import loompy as lp
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad

os.chdir('/data1/ivanir/Ilaria2021/data')

adata = sc.read('structure-identity.loom', sparse=True)

file = open('cell_barcodes_scvelo.txt', 'r')
barcodes = file.read().splitlines()

adata.obs_names_make_unique()
adata.var_names_make_unique()

adata = adata[barcodes]
adata.obs = pd.read_csv('meta_annotated_updated_scvelo.csv', low_memory=False)

scv.pp.pca(adata, n_comps=30)
pca = pd.read_csv('pca_corrected_seurat_updated_scvelo.csv', sep =',', low_memory=False)
pca = pd.DataFrame(pca)
adata.obsm['X_pca'] = pca.values
adata.obsm['X_pca'] = adata.obsm['X_pca'].astype(float)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2500)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.umap(adata)

umap = pd.read_csv('umap_layout_batch_corrected_contamination_cleaned.csv', sep =',', low_memory=False)
umap = pd.DataFrame(umap)
adata.obsm['X_umap'] = umap.values
adata.obsm['X_umap'] = adata.obsm['X_umap'].astype(float)

adata.obs[['scDblFinder.originAmbiguous']] = adata.obs[['scDblFinder.originAmbiguous']].astype('string')

scv.tl.recover_dynamics(adata, n_jobs=60)
scv.tl.velocity(adata, mode='dynamical')
adata.write('./write/scvelo_analysis_dynamical_contamination_cleaned.h5ad')

scv.tl.velocity_graph(adata, n_jobs=60)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype_annotation')

celltype_annotation_colours = {
"Mitotic RG" : "#005dd8",
"Cycling RG" : "#4aadd6",
"Differentiating RG" : "#6083d1",
"Ventricular RG" : "#02d7ec",
"Glycolytic RG" : "#0086cb",
"Transcriptionally active RG" : "#4aac8d",
"Choroid plexus" : "#935de6",
"Cortical hem" : "#9e4d70",
"IPC" : "#ff8ba6",
"Committed neurons" : "#90b600",
"Inhibitory neurons" : "#228B22",
"Cajal Retzius cells" : "#00ff00",
"Migrating excitatory neurons" : "#a6023e",
"UL enriched neurons" : "#fa2274",
"DL enriched neurons" : "#9e5d56",
"Migrating excitatory neurons" : "#a6023e",
"Mature excitatory neurons" : "#d3b000",
}

adata.uns['celltype_annotation_colors'] = [celltype_annotation_colours[i] for i in sorted(np.unique(adata.obs['celltype_annotation']))]

scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype_annotation')

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, basis='umap')

adata.write('./write/scvelo_analysis_contamination_cleaned.h5ad')
adata.obs.to_csv('structure-identity_scvelo_metadata_contamination_cleaned.csv')

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]

scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=100, yticklabels=True)

os.chdir('/data1/ivanir/Ilaria2021/data')

adata_tmp = sc.read('./write/scvelo_analysis_contamination_cleaned.h5ad')

adata = sc.read('structure-identity.loom', sparse=True)

file = open('cell_barcodes_scvelo.txt', 'r')
barcodes = file.read().splitlines()

adata.obs_names_make_unique()
adata.var_names_make_unique()

adata = adata[barcodes]
adata.obs = pd.read_csv('meta_annotated_updated_scvelo.csv', low_memory=False)

scv.pp.pca(adata, n_comps=30)
pca = pd.read_csv('pca_corrected_seurat_updated_scvelo.csv', sep =',', low_memory=False)
pca = pd.DataFrame(pca)
adata.obsm['X_pca'] = pca.values
adata.obsm['X_pca'] = adata.obsm['X_pca'].astype(float)

adata.obs['latent_time'] = adata_tmp.obs['latent_time'].values

celltype_annotation_colours = {
"Mitotic RG" : "#005dd8",
"Cycling RG" : "#4aadd6",
"Differentiating RG" : "#6083d1",
"Ventricular RG" : "#02d7ec",
"Glycolytic RG" : "#0086cb",
"Transcriptionally active RG" : "#4aac8d",
"Choroid plexus" : "#935de6",
"Cortical hem" : "#9e4d70",
"IPC" : "#ff8ba6",
"Committed neurons" : "#90b600",
"Inhibitory neurons" : "#228B22",
"Cajal Retzius cells" : "#00ff00",
"Migrating excitatory neurons" : "#a6023e",
"UL enriched neurons" : "#fa2274",
"DL enriched neurons" : "#9e5d56",
"Migrating excitatory neurons" : "#a6023e",
"Mature excitatory neurons" : "#d3b000",
}

adata.uns['celltype_annotation_colors'] = [celltype_annotation_colours[i] for i in sorted(np.unique(adata.obs['celltype_annotation']))]

file = open('FCD_up_genes.txt', 'r')
fcd_genes = file.read().splitlines()

file = open('CTRL_up_genes.txt', 'r')
ctrl_genes = file.read().splitlines()

file = open('HSM_genes.txt', 'r')
hsm_genes = file.read().splitlines()

file = open('LSM_genes.txt', 'r')
lsm_genes = file.read().splitlines()

scv.pl.heatmap(adata, var_names=fcd_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=100, yticklabels=True, font_scale= .5)

scv.pl.heatmap(adata, var_names=ctrl_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=100, yticklabels=True, font_scale= .5)

scv.pl.heatmap(adata, var_names=hsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=100, yticklabels=True, font_scale= .5)

scv.pl.heatmap(adata, var_names=lsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=100, yticklabels=True, font_scale= .5)

file = open('HSM_neu_genes.txt', 'r')
hsm_genes = file.read().splitlines()

file = open('LSM_neu_genes.txt', 'r')
lsm_genes = file.read().splitlines()

scv.pl.heatmap(adata, var_names=hsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=500, yticklabels=True, font_scale= .5)

scv.pl.heatmap(adata, var_names=lsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=500, yticklabels=True, font_scale= .5)

hsm = scv.pl.heatmap(adata, var_names=hsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=1000, yticklabels=True, show=False)

lsm = scv.pl.heatmap(adata, var_names=lsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=1000, yticklabels=True, show=False)

hsm.data2d.to_csv(r'hsm_ltime.csv', index = True)
lsm.data2d.to_csv(r'lsm_ltime.csv', index = True)

file = open('HSM_prog_genes.txt', 'r')
hsm_genes = file.read().splitlines()

file = open('LSM_prog_genes.txt', 'r')
lsm_genes = file.read().splitlines()

hsm = scv.pl.heatmap(adata, var_names=hsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=1000, yticklabels=True, show=False)

lsm = scv.pl.heatmap(adata, var_names=lsm_genes, sortby='latent_time', col_color='celltype_annotation', n_convolve=1000, yticklabels=True, show=False)

hsm.data2d.to_csv(r'hsm_ltime_prog.csv', index = True)
lsm.data2d.to_csv(r'lsm_ltime_prog.csv', index = True)

adata_neu = adata[adata.obs['celltype_annotation'].isin(['Committed neurons', 'Inhibitory neurons', 'Cajal Retzius cells', 'Migrating excitatory neurons', 'UL enriched neurons', 'DL enriched neurons', 'Migrating excitatory neurons','Mature excitatory neurons'])]
adata_neu.obs[['scDblFinder.originAmbiguous']] = adata_neu.obs[['scDblFinder.originAmbiguous']].astype('string')




