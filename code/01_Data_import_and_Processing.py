import os
import sys
import numpy as np
import pandas as pd

import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt

cellranger_dir="/home/yangyu/Project/skin-aging/cellranger/"
df = pd.read_table('sample_info.txt')

#Load
filenames = [os.path.join(cellranger_dir, i,'outs/filtered_feature_bc_matrix') for i in df['sample_id'].tolist()]

# Concatenate and save
adatas = [sc.read_10x_mtx(filename, cache=True) for filename in filenames]

adata = adatas[0].concatenate(adatas[1:], join='inner', batch_key='sample_id', batch_categories=df['sample_id'].tolist(), index_unique='-')

save_file = "1.adata_merge_raw.h5ad"
adata.write(save_file)

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata.var['ribo'] = adata.var_names.str.startswith('RPS','RPL')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)

adata.write('1.1.adata_merge_raw_filtered.h5ad')

sc.pl.violin(adata, ['n_genes', 'total_counts', 'pct_counts_mt','pct_counts_ribo'],
             jitter=0.4, multi_panel=True,save='_raw.pdf')

adata = adata[adata.obs.n_genes_by_counts < 8000, :]
adata = adata[adata.obs.n_genes_by_counts > 1000, :]

adata = adata[adata.obs.pct_counts_mt < 20, :]
adata
adata.write('1.2.adata_merge_filtered.h5ad')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=False,
    flavor="seurat_v3",
    batch_key="sample_id"
)

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True,save='.pdf')

sc.pp.neighbors(adata, n_neighbors=40, n_pcs=40)
sc.tl.leiden(adata,resolution=0.3,key_added="leiden_res0.3")
sc.tl.umap(adata)

sc.tl.leiden(adata,resolution=0.1,key_added="leiden_res0.1")

sc.pl.umap(adata, color=['leiden_res0.1'],save='_res0.1.pdf')

sc.pl.umap(adata, color=['leiden_res0.3','sample_group','sample_id'],ncols=1,frameon=False,save='_merge_0.3.pdf')

marker_genes_dict = {
    #KC
    'BAS': ["KRT14","KRT5"],
    'SPN': ["KRT10","KRT1"],
    'GRN': ["FLG","KRT2","IVL","KRT20"],
    'Sebocyte':["FADS2","PPARG","FAR2","KRT7","MC5R","KRT4","MUC1"],#"MGST1",
    'HFC' : ["SOX9", "KRT6B", "SFRP1",'LHX2','CD34','TCF3','NFATC1'], #
    'Epidermal':['EPCAM','KRT18','KRT8','KRT7','KRT19'],
    'Melanocytes': ["MLANA","DCT","TYRP1","PMEL","GPR143"],
    'Sweat gland cell': ["AQP5","SCGB2A2","SCGB1D2","CA2","MUCL1"],
    'Fb': ["THY1","PDGFRA","LUM", "DCN","COL1A2"],
    'Nc':["NRXN1","SCN7A","NGFR"],#neuron
    'Schwann': ["S100B","PLP1",'MPZ',"SOX10","SOX2"],#"PMP22",
    'Pc-vSMC': ["ACTA2","RGS5","TAGLN","MYLK"],#Pericyte
    'VEC': ["PECAM1","VWF"],
    'LEC': ["LYVE1","PROX1"],
    'T': ["CD3D","CD3E"],
    'Mast': ["TPSAB1","TPSB2","CMA1"],
    'LC': ["CD1A","CD207"],
    'Macrophage':["CD68","CD163"],
    'Monocyte':["CD14"],
    'B': ["CD19","MS4A1","CD79A", "CD79B", "BLNK","CD5"] #MS4A1
}
sc.pl.umap(adata, color=['leiden_res0.3','cell_type','sample_group','sample_id'],frameon=False,legend_loc='on data',legend_fontweight="normal",
           legend_fontsize='x-small',ncols=2,save='_celltype.pdf')
adata.write('1.2_skin_raw.h5ad')
