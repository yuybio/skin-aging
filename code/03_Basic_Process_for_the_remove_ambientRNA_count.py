import os
import sys
import numpy as np
import pandas as pd

import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

decontx_dir="/home/yangyu/Project/skik-aging/2.ambiantRNAremoved_DecontX/h5adFiles/"
df = pd.read_table('sample_info.txt')
samples = df['sample_id'].tolist()

#Load
filenames = [decontx_dir+i+'.h5ad' for i in samples]

# Concatenate and save
adatas = [sc.read_h5ad(filename) for filename in filenames]

adata = adatas[0].concatenate(adatas[1:], join='inner', batch_key='sample_id', batch_categories=samples, index_unique='-')

save_file = "2.adata_merge_DecontX.h5ad"
adata.write(save_file)

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Doublet_dection_scrublet median -/+ two median absolute deviations of the doublet score

RUNs, DSs, CELLs, THRs, MEDs, MADs, CUTs, no_thr = [], [], [], [], [], [], [], []
# Loop through channels in anndata object:
orig_stdout = sys.stdout
sys.stdout = open('1.scrublet_output_file_mad.txt', 'w')

for run in adata.obs['sample_id'].unique():
    print(run)
    ad = adata[adata.obs['sample_id'] == run, :]
    x = ad.layers['counts'] # use raw count
    scrub = scr.Scrublet(x)
    ds, prd = scrub.scrub_doublets()
    RUNs.append(run)
    DSs.append(ds)
    CELLs.append(ad.obs_names)
    # MAD calculation of threshold:
    MED = np.median(ds)
    MAD = robust.mad(ds)
    CUT = (MED + (4 * MAD))
    MEDs.append(MED)
    MADs.append(MAD)
    CUTs.append(CUT)

    try:  # not always can calculate automatic threshold
        THRs.append(scrub.threshold_)
        print('Threshold found by scrublet')
    except:
        THRs.append(0.4)
        no_thr.append(run)
        print('No threshold found, assigning 0.4 to', run)
        scrub.call_doublets(threshold=0.4) # so that it can make the plot
    fig = scrub.plot_histogram()
    fig[0].savefig(run + '_scrublet_histogram.png')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_))
    fig2=scrub.plot_embedding('UMAP', order_points=True)
    fig2[0].savefig(run + '_scrublet_umap.png')

    # Alternative histogram for MAD-based cutoff
    scrub.call_doublets(threshold=CUT)
    fig = scrub.plot_histogram()
    fig[0].savefig(run + '_mad_' + 'scrublet_histogram.png')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_))
    fig2=scrub.plot_embedding('UMAP', order_points=True)
    fig2[0].savefig(run + '_mad_' +'scrublet_umap.png')
    plt.close('all')
    print()
    print()

print()
print('The following sample(s) do not have automatic threshold:')
print(no_thr)

sys.stdout.close()
sys.stdout = orig_stdout

ns = np.array(list(map(len, DSs)))

tbl = pd.DataFrame({
    'sample': np.repeat(RUNs, ns),
    'ds': np.concatenate(DSs),
    'thr': np.repeat(THRs, ns),
    'mad_MED': np.repeat(MEDs, ns),
    'mad_MAD': np.repeat(MADs, ns),
    'mad_thr': np.repeat(CUTs, ns),
    }, index=np.concatenate(CELLs))

tbl['auto_prd'] = tbl['ds'] > tbl['thr']
tbl['mad_prd'] = tbl['ds'] > tbl['mad_thr']

tbl.to_csv('1.doublets_score_mad.csv', header=True, index=True)

adata.obs['auto_prd']=tbl['auto_prd']
adata.obs['mad_prd']=tbl['mad_prd']


adata.obs['auto_prd'] = adata.obs['auto_prd'].replace(True, 'dobulet')
adata.obs['auto_prd'] = adata.obs['auto_prd'].replace(False, 'singlet')
adata.obs['mad_prd'] = adata.obs['mad_prd'].replace(True, 'dobulet')
adata.obs['mad_prd'] = adata.obs['mad_prd'].replace(False, 'singlet')

sc.pl.violin(adata, ['n_genes', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True,save='_raw.pdf')

adata = adata[adata.obs.n_genes_by_counts < 8000, :]
adata = adata[adata.obs.n_genes_by_counts > 1000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

adata.layers["decontXcounts"] = adata.X.copy()
adata.write(h5ad_dir+'/'+'2.2.adata_ambiantRNAremoved_filtered.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=False,
    layer="decontXcounts",
    flavor="seurat_v3",
    batch_key="sample_id"
)

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True,save='.pdf')

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)
sc.tl.leiden(adata,resolution=0.3,key_added="leiden_res0.3")
sc.tl.umap(adata)

sc.pl.umap(adata, color=['leiden_res0.3','sample_id','sample_group'],ncols=3,frameon=False,legend_loc='on data',legend_fontweight="normal",
           legend_fontsize='x-small',save='_merge_0.3.pdf')

marker_genes_dict = {
    'BAS': ["KRT14","KRT5"],
    'SPN': ["KRT10","KRT1"],
    'GRN': ["FLG","KRT2","IVL","KRT20"],
    'Sebocyte':["FADS2","PPARG","FAR2","KRT7","MC5R","KRT4","MUC1"],#"MGST1",
    'HFC' : ["SOX9", "KRT6B", "SFRP1",'LHX2','CD34','TCF3','NFATC1'],
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

adata.write('2.3_atlas_ambiantRNAremoved.h5ad')
