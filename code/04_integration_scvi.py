import os
import sys
import pandas as pd
import scanpy as sc
import numpy as np
import scvi

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore
from scipy.sparse import spmatrix
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance
from scipy.io import mmread
from scipy.sparse import issparse, csr_matrix
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform

from plotnine import *
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors
import mplscience

adata = sc.read_h5ad('2.2.adata_ambiantRNAremoved_filtered.h5ad')

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=False,
    layer="decontXcounts",
    flavor="seurat_v3",
    batch_key="sample_id"
)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="decontXcounts",
    #batch_key = 'sample_id',
    categorical_covariate_keys=['stage'],
    continuous_covariate_keys=["pct_counts_mt","total_counts"]
)

model = scvi.model.SCVI(adata,
                        gene_likelihood = 'nb',
                        n_latent=50,
                        n_layers=3,
                        dropout_rate=0.1,
                        dispersion = "gene")

%%time
model.train()
model.save("my_model/")

%%time
latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)


sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True,n_pcs=50,save='.pdf')

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI",n_pcs=40, n_neighbors=20)
sc.tl.umap(adata)
sc.tl.leiden(adata,resolution=0.1,key_added="leiden_scVI_res0.1")

sc.tl.leiden(adata,resolution=0.3,key_added="leiden_scVI_res0.3")
sc.tl.leiden(adata,resolution=1,key_added="leiden_scVI_res1")
sc.tl.leiden(adata,resolution=1.5,key_added="leiden_scVI_res1.5")

sc.pl.umap(adata, color=['leiden_scVI_res1.5','leiden_scVI_res1','leiden_scVI_res0.3'], frameon=False,
           legend_loc='on data',legend_fontweight="normal",
           legend_fontsize='x-small')

marker_genes_dict = {
    'BAS': ["KRT14","KRT5"],
    'SPN': ["KRT10","KRT1"],
    'GRN': ["FLG","IVL"],
    'Sebocyte':["FADS2","PPARG","FAR2","KRT7","MC5R","KRT4","MUC1"],
    'HFC' : ["SOX9", "KRT6B", 'LHX2','CD34','TCF3','NFATC1'],
    'Epidermal':['EPCAM','KRT18','KRT8','KRT7','KRT19'],
    'Melanocytes': ["MLANA","DCT","TYRP1","PMEL","GPR143"],
    'Sweat gland cell': ["AQP5","SCGB1D2","CA2","MUCL1"],
    'Fb': ["THY1","PDGFRA","LUM", "DCN","COL1A2"],
    'Nc':["NRXN1","SCN7A","NGFR"],#neuron
    'Sc': ["S100B","PLP1",'MPZ',"SOX10","SOX2"],
    'Pc-vSMC': ["ACTA2","RGS5","TAGLN","MYLK"],
    'VEC': ["PECAM1","VWF"],
    'LEC': ["LYVE1","PROX1"],
    'T': ["CD3D","CD3E"],
    'Mast': ["TPSAB1","TPSB2","CMA1"],
    'LC': ["CD1A","CD207"],
    'Macrophage':["CD68","CD163"],
    'Monocyte':["CD14"],
    'B': ["CD19","MS4A1","CD79A", "CD79B", "BLNK","CD5"]
}

adata
dp=sc.pl.dotplot(adata, marker_genes_dict,'leiden_scVI_res1.5',color_map="RdBu_r", return_fig=True,vmin=-1, vmax=3)
dp.add_totals(color='#757575')
dp.savefig('figures/scVI_res1.5_dotplot.pdf',dpi=1080,bbox_inches='tight')



pseudobulk_count = pd.DataFrame(adata.layers['counts'].A, 
                                        index=adata.obs.index, 
                                        columns=adata.var_names).groupby(adata.obs.sample_id).sum(0)

pseudobulk_count.T.to_csv('0.all_celltype/whole_skin_pseudobulk_rawCount.csv')

for cell_type in adata.obs['cell_type'].unique():
    adata_subset = adata[adata.obs['cell_type'] == cell_type, :].copy()
    pseudobulk_count = pd.DataFrame(
        adata_subset.layers['counts'].A,
        index=adata_subset.obs.index,
        columns=adata_subset.var_names
    ).groupby(adata_subset.obs['sample_id']).sum(0)

    output_path = f'/0.all_celltype/{cell_type}_pseudobulk_rawCount.csv'
    pseudobulk_count.T.to_csv(output_path)  

    print(f"{cell_type} pseudobulk, {output_path}")
