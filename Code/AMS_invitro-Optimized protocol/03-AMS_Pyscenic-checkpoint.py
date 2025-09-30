#!/usr/bin/env python
# coding: utf-8

# # 0. Preprocessing

# In[1]:


import anndata as ad
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# In[2]:


import json
import zlib
import base64
import loompy as lp
import seaborn as sns


# In[3]:


import os

# Set up global parameters
OS_path = "/Users/gclyu07/Desktop/AMS_analysis/"
OS_path_datasets = os.path.join(OS_path, "dataset/")
OS_path_inputs = os.path.join(OS_path, "inputs/")
OS_path_outputs = os.path.join(OS_path, "outputs/")

seed = 1121


# In[15]:


BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"


# In[4]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/AMS_analysis/outputs/s5_0728_AMS_with_unspliced.h5ad")


# In[6]:


adata


# In[5]:


cell_list = {
    0:'Rgl1',
    1:'Rgl2',
    2:'Rgl3',
    3:'ProgFP',
    4:'NProg',
    5:'NbM',
    6:'GabaNb',
    7:'DA'}
adata.obs['bimod_oc'] = (
adata.obs['bimod_oc']
.astype('category')
.map(cell_list)
)


# In[6]:


cell_list = {
    0: 'Rgl1',
    1: 'Rgl2',
    2: 'Rgl3',
    3: 'ProgFP.S',
    4: 'ProgFP.G2M',
    5: 'ProgFP.Development',
    6: 'NProg',
    7: 'NbM',
    8: 'GabaNb',
    9: 'DA.Neuron projection',
    10: 'DA.Synaptic assembly',
}
adata.obs['bimod_oc_des'] = (
adata.obs['bimod_oc_des']
.astype('category')
.map(cell_list)
)


# In[41]:


adata.obsm['umap'] = adata.obsm['X_umap_pca_RNA_0.3']


# In[42]:


sc.set_figure_params(figsize=(4, 3))
sc.pl.umap(adata, color='orig.ident')
sc.pl.umap(adata, color='bimod_oc_des')
sc.pl.umap(adata, color='bimod_oc')


# In[43]:


sc.set_figure_params(figsize=(10, 8))
sc.pl.umap(adata, color = ['orig.ident'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', 
           palette=['#E5D2DD', '#5F3D69', '#C5DEBA'])


# In[44]:


sc.set_figure_params(figsize=(10, 8))
sc.pl.umap(adata, color = ['bimod_oc'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='')


# In[45]:


sc.set_figure_params(figsize=(10, 8))
sc.pl.umap(adata, color = ['bimod_oc_des'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='')


# In[18]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(2, 1.5))
sc.pl.umap(adata, color=["FOXA2", "LMX1A", "EN1", "NR4A2", "TH", "DDC", "OTX2", "SOX6", "LMO3", 
                     "PBX1","KCNJ6","SHH", "ASCL1", "NEUROG2", "DCX","SLC18A2", "bimod_oc_des"],
           title=["FOXA2", "LMX1A", "EN1", "NR4A2", "TH"," DDC", "OTX2", "SOX6", "LMO3", 
                     "PBX1","KCNJ6","SHH", "ASCL1", "NEUROG2", "DCX","SLC18A2", "Cell types"], ncols=3, use_raw=False, cmap='magma_r')


# In[7]:


# collect SCENIC AUCell output
lf = lp.connect( "/Users/gclyu07/Desktop/AMS_analysis/outputs/AMS_scenic.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


# # 1. Regulon projections

# In[8]:


ad_auc_mtx = ad.AnnData(auc_mtx)


# In[16]:


sc.pp.neighbors(ad_auc_mtx, n_neighbors=10, metric="correlation")
sc.tl.umap(ad_auc_mtx)
sc.tl.tsne(ad_auc_mtx)


# In[17]:


adata.obsm["X_umap_aucell"] = ad_auc_mtx.obsm["X_umap"]
adata.obsm["X_tsne_aucell"] = ad_auc_mtx.obsm["X_tsne"]


# In[25]:


sc.set_figure_params(figsize=(6, 4))
sc.pl.embedding(adata, basis="X_umap_aucell", color="orig.ident", palette=['#E5D2DD', '#5F3D69', '#C5DEBA'],
                add_outline=True, legend_fontsize=10, legend_fontoutline=2,frameon=False, size=25, title = "AUCell UMAP")


# In[26]:


sc.set_figure_params(figsize=(6, 4))
sc.pl.embedding(adata, basis="X_umap_aucell", color="bimod_oc", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell UMAP")


# In[27]:


sc.set_figure_params(figsize=(6, 4))
sc.pl.embedding(adata, basis="X_umap_aucell", color="bimod_oc_des", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell UMAP")


# # 3. Regulon clustering

# In[28]:


auc_mtx["bimod_oc_des"] = adata.obs['bimod_oc_des'] 
mean_auc_by_cell_type = auc_mtx.groupby("bimod_oc_des").mean()


# In[29]:


top_n = 50
top_tfs = mean_auc_by_cell_type.max(axis=0).sort_values(ascending=False).head(top_n)
mean_auc_by_cell_type_top_n = mean_auc_by_cell_type[
    [c for c in mean_auc_by_cell_type.columns if c in top_tfs]
]


# In[30]:


tf_names = top_tfs.index.str.replace(r"\(\+\)", "", regex=True)
adata_batch_top_tfs = adata[:, adata.var_names.isin(tf_names)]


# In[31]:


sc.pl.matrixplot(
    adata,
    tf_names,
    use_raw=False,
    groupby="bimod_oc_des",
    cmap="Reds",
    dendrogram=True,
    figsize=[15, 4],
    standard_scale="group",
)


# In[32]:


scv.set_figure_params('scvelo') 
sns.clustermap(
    mean_auc_by_cell_type_top_n,
    figsize=[12, 6],
    cmap="viridis_r",
    xticklabels=True,
    yticklabels=True,
    standard_scale=0,
    cbar_pos=(1, .3, .01, .4),
)


# In[33]:


scv.set_figure_params('scvelo') 
sns.clustermap(
    mean_auc_by_cell_type_top_n,
    figsize=[12, 6],
    cmap="viridis_r",
    xticklabels=True,
    yticklabels=True,
    standard_scale=1,
    cbar_pos=(1, .3, .01, .4),
)


# ______

# In[36]:


auc_mtx["bimod_oc"] = adata.obs['bimod_oc'] 
mean_auc_by_cell_type = auc_mtx.groupby("bimod_oc").mean()


# In[37]:


top_n = 50
top_tfs = mean_auc_by_cell_type.max(axis=0).sort_values(ascending=False).head(top_n)
mean_auc_by_cell_type_top_n = mean_auc_by_cell_type[
    [c for c in mean_auc_by_cell_type.columns if c in top_tfs]
]


# In[38]:


tf_names = top_tfs.index.str.replace(r"\(\+\)", "", regex=True)
adata_batch_top_tfs = adata[:, adata.var_names.isin(tf_names)]


# In[39]:


sc.pl.matrixplot(
    adata,
    tf_names,
    use_raw=False,
    groupby="bimod_oc_des",
    cmap="Reds",
    dendrogram=True,
    figsize=[15, 4],
    standard_scale="group",
)


# In[43]:


scv.set_figure_params('scvelo') 
sns.clustermap(
    mean_auc_by_cell_type_top_n,
    figsize=[12, 5],
    cmap="viridis_r",
    xticklabels=True,
    yticklabels=True,
    standard_scale=0,
    cbar_pos=(1, .3, .01, .4),
)


# In[42]:


scv.set_figure_params('scvelo') 
sns.clustermap(
    mean_auc_by_cell_type_top_n,
    figsize=[12, 5],
    cmap="viridis_r",
    xticklabels=True,
    yticklabels=True,
    standard_scale=1,
    cbar_pos=(1, .3, .01, .4),
)


# In[ ]:





# # 3. Visualization from pyscenic

# In[9]:


from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss


# In[10]:


def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(pal),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f


# In[11]:


# collect SCENIC AUCell output
lf = lp.connect( "/Users/gclyu07/Desktop/AMS_analysis/outputs/AMS_scenic.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


# In[48]:


bin_mtx, thresholds = binarize(auc_mtx,seed=1121,num_workers=10)

fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(8, 5), dpi=100)

plot_binarization(auc_mtx, 'FOXA2(+)', thresholds['FOXA2(+)'], ax=ax1)
plot_binarization(auc_mtx, 'EN1(+)', thresholds['EN1(+)'], ax=ax2)
plot_binarization(auc_mtx, 'OTX2(+)', thresholds['OTX2(+)'], ax=ax3)
plot_binarization(auc_mtx, 'SOX6(+)', thresholds['SOX6(+)'], ax=ax4)
plot_binarization(auc_mtx, 'NEUROD1(+)', thresholds['NEUROD1(+)'], ax=ax5)
plot_binarization(auc_mtx, 'NHLH1(+)', thresholds['NHLH1(+)'], ax=ax6)

# In[50]:


from matplotlib.colors import to_hex
N_COLORS = len(adata.obs['bimod_oc_des'].dtype.categories)
cmap = plt.get_cmap('tab20')
COLORS = [to_hex(cmap(i / N_COLORS)) for i in range(N_COLORS)]
cell_type_color_lut = dict(zip(adata.obs['bimod_oc_des'].dtype.categories, COLORS))
#cell_type_color_lut = dict(zip(adata.obs.cell_type.dtype.categories, adata.uns['cell_type_colors']))
cell_id2cell_type_lut = adata.obs['bimod_oc_des'].to_dict()

bw_palette = sns.xkcd_palette(["white", "black"])


# In[51]:


adata.obs['bimod_oc_des'].dtype.categories


# In[52]:


adata.obs['bimod_oc_des'].to_dict()


# In[53]:


sns.set()
sns.set_style("whitegrid")
fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])


# In[54]:


sns.set()
sns.set(font_scale=0.7)
fig = palplot(sns.color_palette(COLORS), adata.obs['bimod_oc_des'].dtype.categories, size=2)


# In[55]:


sns.set()
sns.set(font_scale=1.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
g = sns.clustermap(bin_mtx.T, 
               col_colors=auc_mtx.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
               cmap=bw_palette, figsize=(20,20))
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cell Type'])
g.cax.set_visible(False)


# ## 3.2 Generate sequence logos

# In[12]:


import operator as op
import cytoolz


# In[16]:


def fetch_logo(regulon, base_url = BASE_URL):
    for elem in regulon.context:
        if elem.endswith('.png'):
            return '<img src="{}{}" style="max-height:124px;"></img>'.format(base_url, elem)
    return ""


# In[17]:


def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)


# In[18]:


lf = lp.connect( "/Users/gclyu07/Desktop/AMS_analysis/outputs/AMS_scenic.loom", mode='r+', validate=False )


# In[19]:


# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  r[r==1].index.values


# In[20]:


def filter_regulons(motifs, db_names=("hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings",)):

    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(cytoolz.compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0) 
                                                                      & ((motifs['Annotation'] == 'gene is directly annotated')
                                                                        | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                                           & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))


# In[21]:


df_motifs = load_motifs("outputs/AMS_reg.csv")
regulons = filter_regulons(df_motifs)


# In[65]:


df_regulons = pd.DataFrame(data=[map(op.attrgetter('name'), regulons),
                                 map(len, regulons),
                                 map(fetch_logo, regulons)], index=['name', 'count', 'logo']).T


# In[64]:


from IPython.display import HTML, display


# In[66]:


MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
display(HTML(df_regulons.head().to_html(escape=False)))
pd.set_option('display.max_colwidth', MAX_COL_WIDTH)


# ## 3.3 CELL TYPE SPECIFIC REGULATORS - Z-SCORE

# In[22]:


add_scenic_metadata(adata, auc_mtx, regulons)


# In[68]:


df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['bimod_oc_des']]
df_results = ((df_scores.groupby(by='bimod_oc_des').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results[(df_results.Z >= 1.0)].sort_values('Z', ascending=False)


# In[69]:


df_results['regulon'] = df_results['regulon'].str.replace(r'\(\+\)', '', regex=True)


# In[70]:


df_results


# In[71]:


df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 1.1].sort_values('Z', ascending=False),
                           index='bimod_oc_des', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(25, 15))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=0.5, cbar=False, square=True, linecolor='grey', 
            cmap="YlGnBu", annot_kws={"size": 8},)
ax1.set_ylabel('')
ax1.set_xlabel('')


# In[73]:


df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['bimod_oc']]
df_results = ((df_scores.groupby(by='bimod_oc').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results[(df_results.Z >= 1.0)].sort_values('Z', ascending=False)


# In[74]:


df_results['regulon'] = df_results['regulon'].str.replace(r'\(\+\)', '', regex=True)


# In[75]:


df_results


# In[77]:


df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 1.1].sort_values('Z', ascending=False),
                           index='bimod_oc', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(25, 15))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=0.5, cbar=False, square=True, linecolor='grey', 
            cmap="YlGnBu", annot_kws={"size": 8},)
ax1.set_ylabel('')
ax1.set_xlabel('')


# ## 3.4 CELL TYPE SPECIFIC REGULATORS - RSS

# In[78]:


# collect SCENIC AUCell output
lf = lp.connect( "/Users/gclyu07/Desktop/AMS_analysis/outputs/AMS_scenic.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


# In[79]:


rss = regulon_specificity_scores(auc_mtx, adata.obs.bimod_oc_des)
rss.head()


# In[80]:


adata.obs.bimod_oc_des.unique()


# In[84]:


sns.set()
sns.set(style='whitegrid', font_scale=0.8)
fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(9, 3), dpi=400)
plot_rss(rss, 'DA.Neuron projection', ax=ax1)
ax1.set_xlabel('')
plot_rss(rss, 'DA.Synaptic assembly', ax=ax2)
ax2.set_xlabel('')
ax2.set_ylabel('')
plot_rss(rss, 'NbM', ax=ax3)
ax5.set_xlabel('')
ax5.set_ylabel('')


# In[85]:


rss = regulon_specificity_scores(auc_mtx, adata.obs.bimod_oc)
rss.head()


# In[86]:


sns.set()
sns.set(style='whitegrid', font_scale=0.8)
fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(1, 4, figsize=(12,5), dpi=400)
plot_rss(rss, 'DA', ax=ax1)
ax1.set_xlabel('')
plot_rss(rss, 'NbM', ax=ax2)
ax2.set_xlabel('')
ax2.set_ylabel('')
plot_rss(rss, 'ProgFP', ax=ax3)
ax3.set_xlabel('')
ax3.set_ylabel('')
plot_rss(rss, 'Rgl1', ax=ax4)
ax4.set_xlabel('')
ax4.set_ylabel('')


# In[54]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 2))
sc.pl.umap(adata, color=['UNCX', 'POU2F2', 'PBX1', 'PITX3','ASCL1','bimod_oc_des'],
           title=['UNCX', 'POU2F2', 'PBX1', 'PITX3','ASCL1','Cell types'], ncols=3, use_raw=False, cmap='magma_r')


# In[51]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 2))
sc.pl.umap(adata, color=['Regulon(UNCX(+))', 'Regulon(POU2F2(+))', 'Regulon(PBX1(+))',
                         'Regulon(PITX3(+))','Regulon(ASCL1(+))','bimod_oc_des'],
           title=['UNCX', 'POU2F2', 'PBX1', 'PITX3','ASCL1','Cell types'], ncols=3, use_raw=False, cmap='viridis_r')


# ## 3.5 AUCELL data

# In[89]:


aucell_adata = sc.AnnData(X=auc_mtx.sort_index(),dtype=np.float32)


# In[90]:


aucell_adata.obs = adata.obs


# In[91]:


import re
regulon_columns = [col for col in aucell_adata.obs.columns if col.startswith("Regulon")]
new_column_names = [re.search(r'Regulon\((.*?)\(\+\)', col).group(1) for col in regulon_columns]


# In[92]:


rename_mapping = dict(zip(regulon_columns, new_column_names))


# In[93]:


aucell_adata.obs.rename(columns=rename_mapping, inplace=True)


# In[94]:


df = df_results[(df_results.Z >= 1.1)].sort_values('Z', ascending=False)
names = df.iloc[:, 1].unique().tolist()


# In[95]:


names


# In[96]:


sc.pl.stacked_violin(aucell_adata, names, groupby='bimod_oc_des')


# # 4. SCENIC Visualization 

# In[11]:


library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(dplyr)
library(Seurat)
library(cowplot)


# In[2]:


loom <- open_loom('~/Desktop/AMS_analysis/outputs/AMS_scenic.loom')


# In[3]:


regulons_incidMat <- get_regulons(loom, column.attr.name = 'Regulons')


# In[4]:


regulonAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')


# In[5]:


sec <- readRDS('~/Desktop/AMS_analysis/outputs/s5_0728_AMS_with_unspliced.rds')


# In[6]:


cellinfo <- sec@meta.data[,c('bimod_oc','bimod_oc_des','orig.ident')]
colnames(cellinfo) <- c('bimod_oc','bimod_oc_des','orig.ident')


# ## 4.1 RSS in bimod_oc_des 

# In[7]:


cellTypes <-  as.data.frame(subset(cellinfo,select = 'bimod_oc_des'))
selectedResolution <- "bimod_oc_des"
cellAnnotation = cellTypes[colnames(regulonAUC),
                           selectedResolution]
cellAnnotation = na.omit(cellAnnotation)


# In[8]:


rss <- calcRSS(AUC = getAUC(regulonAUC),
               cellAnnotation = cellAnnotation)
rss = na.omit(rss)


# In[9]:


rssPlot <- plotRSS(rss,
        zThreshold = 2.5,
        cluster_columns = FALSE,
        order_rows = TRUE,
        thr=0.1,
        varName = "cellType",
        col.low = '#330066',
        col.mid = '#66CC66',
        col.high = '#FFCC33')

rssPlot$rowOrder
plotly::ggplotly(rssPlot$plot, width=500, height=1600)


# In[31]:


rss


# In[32]:


Idents(sec) <- 'bimod_oc_des'


# In[35]:


options(repr.plot.width=20, repr.plot.height=10)

plot_list <- list()

for (i in seq(1, 11)) {
  setName <- table(sec@active.ident)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 2) 


# In[12]:


options(repr.plot.width=20, repr.plot.height=5)
VlnPlot(sec, features = c('RFX2','RFX3','ETV6','SOX9'), split.by = 'bimod_oc_des', layer = 'data', ncol=4) + theme(legend.position = "right")

rss <- rss[!rownames(rss) %in% c("BRCA1(+)","MYC(+)","HOXB5(+)","THRB(+)"),]options(repr.plot.width=10, repr.plot.height=20)

plot_list <- list()

for (i in seq(1, 17)) {
  setName <- table(sec$bimod_oc_des)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 5) 
# In[10]:


options(repr.plot.width=4, repr.plot.height=4)

plot_list <- list()
for (i in seq(1, 11)) {
  setName <- table(sec$bimod_oc_des)[i]%>%names()
  print(plotRSS_oneSet(rss[, ], setName = setName, n = 5))
}


# ## 4.2 RSS in bimod_oc

# In[18]:


cellTypes <-  as.data.frame(subset(cellinfo,select = 'bimod_oc'))
selectedResolution <- "bimod_oc"
cellAnnotation = cellTypes[colnames(regulonAUC),
                           selectedResolution]
cellAnnotation = na.omit(cellAnnotation)


# In[19]:


rss <- calcRSS(AUC = getAUC(regulonAUC),
               cellAnnotation = cellAnnotation)
rss = na.omit(rss)


# In[20]:


rssPlot <- plotRSS(rss,
        zThreshold = 2.5,
        cluster_columns = FALSE,
        order_rows = TRUE,
        thr=0.1,
        varName = "cellType",
        col.low = '#330066',
        col.mid = '#66CC66',
        col.high = '#FFCC33')

rssPlot$rowOrder
plotly::ggplotly(rssPlot$plot, width=500, height=1600)


# In[21]:


rss


# In[22]:


Idents(sec) <- 'bimod_oc'


# In[25]:


options(repr.plot.width=10, repr.plot.height=10)

plot_list <- list()

for (i in seq(1, 8)) {
  setName <- table(sec@active.ident)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 3) 

rss <- rss[!rownames(rss) %in% c("BRCA1(+)","MYC(+)","HOXB5(+)","THRB(+)"),]options(repr.plot.width=10, repr.plot.height=10)

plot_list <- list()

for (i in seq(1, 9)) {
  setName <- table(sec$bimod_oc)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 3)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 3) 
# In[27]:


options(repr.plot.width=4, repr.plot.height=5)

plot_list <- list()
for (i in seq(1, 8)) {
  setName <- table(sec$bimod_oc)[i]%>%names()
  print(plotRSS_oneSet(rss[, ], setName = setName, n = 5))
}

