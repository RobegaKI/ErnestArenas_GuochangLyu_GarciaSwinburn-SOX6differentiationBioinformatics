#!/usr/bin/env python
# coding: utf-8

# ###### XMAS_Single Cell Multiome ATAC + Gene Expression_06

# # 0. SCIENIC RUN

# In[ ]:




#!/bin/bash -l

#SBATCH -A naiss2023-22-456
#SBATCH -t 1-00:00:00
#SBATCH -J 051824_XMAS_SCIENC_full
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128gb
#SBATCH --requeue
#SBATCH --output=051824_XMAS_SCIENC_full.out
#SBATCH --error=051824_XMAS_SCIENC_full.err 
#SBATCH --mail-type=END
#SBATCH --mail-user=guochang.lyu@ki.se


## Record the start time
start=`date +%s`
## Record the host being run on
echo "Hostname: $(eval hostname)"

cd /proj/gclyu07d/23_XMAS_SEQ/SCENIC

  call="pyscenic grn ./06_xma_relabeled_monocle3.loom allTFs_hg38.txt -o XMAS_adj.csv --num_workers 20"

## Echo the call
echo $call
## Evaluate the call
eval $call

done


cd /proj/gclyu07d/23_XMAS_SEQ/SCENIC

  call="pyscenic ctx XMAS_adj.csv \
/proj/gclyu07d/23_XMAS_SEQ/SCENIC/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
/proj/gclyu07d/23_XMAS_SEQ/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ./06_xma_relabeled_monocle3.loom \
--output XMAS_reg.csv \
--mask_dropouts \
--num_workers 20"

## Echo the call
echo $call
## Evaluate the call
eval $call

done



cd /proj/gclyu07d/23_XMAS_SEQ/SCENIC

  call="pyscenic aucell ./06_xma_relabeled_monocle3.loom XMAS_reg.csv --output ./07_xmas_scenic.loom --num_workers 20"

## Echo the call
echo $call
## Evaluate the call
eval $call

done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime

~                                                                                                                                    
~                    
# # 1. Preprocessing 

# In[2]:


import os
import anndata
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt


# In[5]:


# test if dict is callable
tel = {'jack': 4098, 'sape': 4139}
list(tel)


# In[3]:


BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"


# In[4]:


adata = sc.read_h5ad('/Users/gclyu07/Desktop/XMAS_analysis/outputs/s03_XMAS_filtered_velocyto_RNA_BiModOC_labeled_with_unspliced.h5ad')


# In[5]:


cell_list = {
0:'Rgl1',
1:'Rgl2',
2:'Rgl3 ECM organization',
3:'Rgl3 Purine metabolism',
4:'Rgl3 Fate commitment',
5:'ProgFP S',
6:'ProgFP G2M',
7:'ProgFP Development',
8:'NbM',
9:'GabaNb',
10:'Gaba',
11:'DA Neuron projection',
12:'DA Synaptic modulation',
13:'DA Synaptic assembly',
14:'DA Neurotransmitter release'}
adata.obs['bimod_oc_des'] = (
adata.obs['bimod_oc_des']
.astype('category')
.map(cell_list)
)


# In[6]:


cell_list = {
    0:'DA',
    1:'NbM',
    2:'Gaba',
    3:'GabaNb',
    4:'ProgFP',
    5:'Rgl1',
    6:'Rgl2',
    7:'Rgl3',}
adata.obs['bimod_oc'] = (
adata.obs['bimod_oc']
.astype('category')
.map(cell_list)
)


# In[8]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['UNCX', 'POU2F2', 'PBX1', 'NEUROD1','NHLH1','bimod_oc_des'],
           title=['UNCX', 'POU2F2', 'PBX1', 'NEUROD1','NHLH1','Cell types'], ncols=3, use_raw=False,
           cmap='magma_r')


# # 2. General visualization 

# In[7]:


import json
import zlib
import base64
import loompy as lp


# In[8]:


# collect SCENIC AUCell output
lf = lp.connect( "/Users/gclyu07/Desktop/XMAS_analysis/outputs/07_XMAS_scenic.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


# In[10]:


ad_auc_mtx = anndata.AnnData(auc_mtx)
sc.pp.neighbors(ad_auc_mtx, n_neighbors=10, metric="correlation")
sc.tl.umap(ad_auc_mtx)
sc.tl.tsne(ad_auc_mtx)


# In[11]:


adata.obsm["X_umap_aucell"] = ad_auc_mtx.obsm["X_umap"]
adata.obsm["X_tsne_aucell"] = ad_auc_mtx.obsm["X_tsne"]


# In[12]:


sc.set_figure_params(figsize=(10, 8))
sc.pl.embedding(adata, basis="X_umap_aucell", color="bimod_oc_des", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell UMAP")


# In[13]:


sc.set_figure_params(figsize=(6, 4))
sc.pl.embedding(adata, basis="X_umap_aucell", color="bimod_oc_des", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell UMAP")


# In[14]:


sc.set_figure_params(figsize=(6, 4))
sc.pl.embedding(adata, basis="X_umap_aucell", color="bimod_oc", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell UMAP")


# In[15]:


sc.set_figure_params(figsize=(10, 8))
sc.pl.embedding(adata, basis="X_tsne_aucell", color="bimod_oc_des", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell tsne")


# In[16]:


sc.set_figure_params(figsize=(6, 4))
sc.pl.embedding(adata, basis="X_tsne_aucell", color="bimod_oc_des", palette='tab20',  add_outline=True, 
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25, title = "AUCell t-sne")


# In[ ]:





# In[18]:


auc_mtx["bimod_oc_des"] = adata.obs["bimod_oc_des"]
mean_auc_by_cell_type = auc_mtx.groupby("bimod_oc_des").mean()


# In[19]:


top_n = 50
top_tfs = mean_auc_by_cell_type.max(axis=0).sort_values(ascending=False).head(top_n)
mean_auc_by_cell_type_top_n = mean_auc_by_cell_type[
    [c for c in mean_auc_by_cell_type.columns if c in top_tfs]
]


# In[20]:


scv.set_figure_params('scvelo') 
sns.clustermap(
    mean_auc_by_cell_type_top_n,
    figsize=[16, 6],
    cmap="viridis_r",
    xticklabels=True,
    yticklabels=True,
    cbar_pos=(1, .3, .01, .4),
)


# In[21]:


tf_names = top_tfs.index.str.replace(r"\(\+\)", "", regex=True)
adata_batch_top_tfs = adata[:, adata.var_names.isin(tf_names)]


# In[22]:


top_tfs.index


# In[23]:


sc.pl.matrixplot(
    adata,
    tf_names,
    groupby="bimod_oc_des",
    cmap="Reds",
    dendrogram=True,
    figsize=[15, 4],
    standard_scale="group",
)


# # 3. pyscenic Visualization 

# ### 3.1 Create heatmap with binarized regulon activity

# In[9]:


from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss


# In[17]:


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


# In[10]:


# collect SCENIC AUCell output
lf = lp.connect( "/Users/gclyu07/Desktop/XMAS_analysis/outputs/07_XMAS_scenic.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


# In[34]:


fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(8, 4), dpi=100)

plot_binarization(auc_mtx, 'FOXA2(+)', thresholds['FOXA2(+)'], ax=ax1)
plot_binarization(auc_mtx, 'EN1(+)', thresholds['EN1(+)'], ax=ax2)
plot_binarization(auc_mtx, 'OTX2(+)', thresholds['OTX2(+)'], ax=ax3)
plot_binarization(auc_mtx, 'SOX6(+)', thresholds['SOX6(+)'], ax=ax4)
plot_binarization(auc_mtx, 'NEUROD1(+)', thresholds['NEUROD1(+)'], ax=ax5)
plot_binarization(auc_mtx, 'NHLH1(+)', thresholds['NHLH1(+)'], ax=ax6)


# In[27]:


bin_mtx, thresholds = binarize(auc_mtx,seed=1121,num_workers=10)


# In[52]:


from matplotlib.colors import to_hex
N_COLORS = len(adata.obs['bimod_oc_des'].dtype.categories)
cmap = plt.get_cmap('tab20')
COLORS = [to_hex(cmap(i / N_COLORS)) for i in range(N_COLORS)]
cell_type_color_lut = dict(zip(adata.obs['bimod_oc_des'].dtype.categories, COLORS))
#cell_type_color_lut = dict(zip(adata.obs.cell_type.dtype.categories, adata.uns['cell_type_colors']))
cell_id2cell_type_lut = adata.obs['bimod_oc_des'].to_dict()

bw_palette = sns.xkcd_palette(["white", "black"])


# In[36]:


adata.obs['bimod_oc_des'].dtype.categories


# In[37]:


adata.obs['bimod_oc_des'].to_dict()


# In[58]:


sns.set()
sns.set_style("whitegrid")
fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])


# In[62]:


sns.set()
sns.set(font_scale=0.7)
fig = palplot(sns.color_palette(COLORS), adata.obs['bimod_oc_des'].dtype.categories, size=2)


# In[63]:


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

# In[11]:


import operator as op
import cytoolz


# In[17]:


def fetch_logo(regulon, base_url = BASE_URL):
    for elem in regulon.context:
        if elem.endswith('.png'):
            return '<img src="{}{}" style="max-height:124px;"></img>'.format(base_url, elem)
    return ""


# In[18]:


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


# In[12]:


lf = lp.connect( "/Users/gclyu07/Desktop/XMAS_analysis/SCENIC/07_XMAS_scenic.loom", mode='r+', validate=False )


# In[13]:


# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  r[r==1].index.values


# In[14]:


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


# In[15]:


df_motifs = load_motifs("outputs/XMAS_reg.csv")
regulons = filter_regulons(df_motifs)

import pickle

# Saving the list to a file
with open('./outputs/regulons.pkl', 'wb') as file:
    pickle.dump(regulons, file)
# In[16]:


df_regulons = pd.DataFrame(data=[map(op.attrgetter('name'), regulons),
                                 map(len, regulons),
                                 map(fetch_logo, regulons)], index=['name', 'count', 'logo']).T


# In[24]:


from IPython.display import HTML, display


# In[25]:


MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
display(HTML(df_regulons.head().to_html(escape=False)))
pd.set_option('display.max_colwidth', MAX_COL_WIDTH)


# ## 3.3 CELL TYPE SPECIFIC REGULATORS - Z-SCORE
import pickle

# Loading the list from the file
with open('./outputs/regulons.pkl', 'rb') as file:
    regulons = pickle.load(file)
# In[17]:


add_scenic_metadata(adata, auc_mtx, regulons)


# In[18]:


df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['bimod_oc_des']]
df_results = ((df_scores.groupby(by='bimod_oc_des').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results[(df_results.Z >= 1.0)].sort_values('Z', ascending=False)


# In[19]:


df_results['regulon'] = df_results['regulon'].str.replace(r'\(\+\)', '', regex=True)


# In[20]:


df_results


# In[21]:


df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 1.1].sort_values('Z', ascending=False),
                           index='bimod_oc_des', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(25, 15))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=0.5, cbar=False, square=True, linecolor='grey', 
            cmap="YlGnBu", annot_kws={"size": 8},)
ax1.set_ylabel('')
ax1.set_xlabel('')


# ## 3.4 CELL TYPE SPECIFIC REGULATORS - RSS

# In[22]:


rss = regulon_specificity_scores(auc_mtx, adata.obs.bimod_oc_des)
rss.head()


# In[53]:


adata.obs.bimod_oc_des.unique()


# In[56]:


sns.set()
sns.set(style='whitegrid', font_scale=0.8)
fig, ((ax1, ax2, ax3, ax4, ax5)) = plt.subplots(1, 5, figsize=(12, 3), dpi=400)
plot_rss(rss, 'DA Neuron projection', ax=ax1)
ax1.set_xlabel('')
plot_rss(rss, 'DA Synaptic assembly', ax=ax2)
ax2.set_xlabel('')
ax2.set_ylabel('')
plot_rss(rss, 'DA Neurotransmitter release', ax=ax3)
ax3.set_xlabel('')
ax3.set_ylabel('')
plot_rss(rss, 'DA Synaptic modulation', ax=ax4)
ax4.set_xlabel('')
ax4.set_ylabel('')
plot_rss(rss, 'NbM', ax=ax5)
ax5.set_xlabel('')
ax5.set_ylabel('')


# In[34]:


rss = regulon_specificity_scores(auc_mtx, adata.obs.bimod_oc)
rss.head()


# In[35]:


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


# In[61]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['UNCX', 'POU2F2', 'PBX1', 'NEUROD1','NHLH1','bimod_oc_des'],
           title=['UNCX', 'POU2F2', 'PBX1', 'NEUROD1','NHLH1','Cell types'], ncols=3, use_raw=False, cmap='magma_r')


# ## 3.5 AUCELL data

# In[23]:


aucell_adata = sc.AnnData(X=auc_mtx.sort_index(),dtype=np.float32)


# In[24]:


aucell_adata.obs = adata.obs


# In[25]:


import re
regulon_columns = [col for col in aucell_adata.obs.columns if col.startswith("Regulon")]
new_column_names = [re.search(r'Regulon\((.*?)\(\+\)', col).group(1) for col in regulon_columns]


# In[26]:


rename_mapping = dict(zip(regulon_columns, new_column_names))
aucell_adata.obs.rename(columns=rename_mapping, inplace=True)


# In[27]:


names = ['AHDC1', 'ARID3A', 'ATF2', 'ATF3', 'ATF4', 'BACH2', 'BARHL2', 'BCL11A', 'BCL3', 'BCLAF1', 'BRCA1', 'CARF', 'CDX2', 'CEBPA', 'CEBPD', 'CERS4', 'CHD1', 'CHD2', 'CRX', 'CUX1', 'DBP', 'DBX2', 'DDIT3', 'DLX2', 'E2F1', 'E2F2', 'E2F3', 'E2F4', 'E2F7', 'E2F8', 'EGR1', 'ELF1', 'ELF2', 'ELF3', 'ELF4', 'ELF5', 'ELK4', 'EMX2', 'EN1', 'EP300', 'ERF', 'ETS1', 'ETV2', 'ETV3', 'ETV6', 'ETV7', 'FEV', 'FOS', 'FOSB', 'FOSL2', 'FOXA1', 'FOXA2', 'FOXJ1', 'FOXK1', 'FOXM1', 'FOXN4', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXP4', 'FOXQ1', 'FOXR1', 'FUBP1', 'GABPB1', 'GATA3', 'GATA4', 'GATA6', 'GCM1', 'GMEB1', 'GSX2', 'HES7', 'HMGA1', 'HMGA2', 'HNF1A', 'HOXA9', 'HOXB5', 'HOXB6', 'HOXC4', 'HOXC8', 'HOXD3', 'IKZF1', 'IKZF2', 'ILF2', 'IRF1', 'IRF4', 'IRF5', 'IRF6', 'IRF7', 'JUN', 'JUNB', 'JUND', 'KLF12', 'KLF17', 'KLF5', 'KLF7', 'KMT2B', 'LCOR', 'LEF1', 'LHX1', 'LHX9', 'MAFA', 'MAFB', 'MAFK', 'MAX', 'MAZ', 'MEOX1', 'MITF', 'MNT', 'MNX1', 'MSX1', 'MTF2', 'MXD1', 'MYBL1', 'MYBL2', 'MYC', 'MYF6', 'MZF1', 'NEUROD1', 'NFATC1', 'NFE2L2', 'NFE2L3', 'NFIA', 'NFIB', 'NFIC', 'NFIX', 'NFYA', 'NFYC', 'NHLH1', 'NKX2-2', 'NKX3-1', 'NR1H4', 'NR2F2', 'NR3C1', 'OTX2', 'OVOL1', 'PAX3', 'PBX1', 'PBX2', 'PDLIM5', 'PITX2', 'POU2F1', 'POU2F2', 'POU3F2', 'POU6F1', 'PRDM1', 'PRDM5', 'RARA', 'REST', 'RFX1', 'RFX2', 'RFX3', 'RFX4', 'RFX5', 'RFXANK', 'RFXAP', 'RHOXF2', 'RUNX3', 'RXRA', 'RXRB', 'SHOX2', 'SMAD1', 'SMAD3', 'SMAD4', 'SND1', 'SOX11', 'SOX12', 'SOX4', 'SOX5', 'SOX6', 'SOX7', 'SP1', 'SP3', 'SPDEF', 'SREBF1', 'SREBF2', 'STAT5A', 'TAF1', 'TAF6', 'TCF12', 'TCF3', 'TCF7L1', 'TCF7L2', 'TEAD4', 'TFAP2B', 'TFDP2', 'TFE3', 'TFEB', 'THRA', 'THRB', 'TP53', 'UNCX', 'USF2', 'VAX1', 'VAX2', 'VDR', 'WT1', 'XBP1', 'YBX1', 'YY1', 'ZBTB4', 'ZBTB7A', 'ZFP14', 'ZFP62', 'ZFX', 'ZFY', 'ZIC1', 'ZNF131', 'ZNF148', 'ZNF235', 'ZNF260', 'ZNF285', 'ZNF343', 'ZNF354C', 'ZNF383', 'ZNF407', 'ZNF432', 'ZNF433', 'ZNF485', 'ZNF518A', 'ZNF519', 'ZNF546', 'ZNF549', 'ZNF582', 'ZNF615', 'ZNF616', 'ZNF626', 'ZNF628', 'ZNF69', 'ZNF692', 'ZNF700', 'ZNF707', 'ZNF711', 'ZNF75A', 'ZNF764', 'ZNF93', 'ZSCAN4']


# In[28]:


sc.pl.stacked_violin(aucell_adata, names, groupby='bimod_oc_des')


# In[29]:


df = df_results[(df_results.Z >= 1.1)].sort_values('Z', ascending=False)
names = df.iloc[:, 1].unique().tolist()


# In[30]:


names


# In[31]:


sc.pl.stacked_violin(aucell_adata, names, groupby='bimod_oc_des')


# # 4. SCENIC Visualization 

# In[19]:


library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(dplyr)
library(Seurat)
library(cowplot)


# In[2]:


loom <- open_loom('~/Desktop/XMAS_analysis/outputs/07_XMAS_scenic.loom')


# In[3]:


regulons_incidMat <- get_regulons(loom, column.attr.name = 'Regulons')


# In[4]:


regulonAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')


# In[5]:


sec <- readRDS('~/Desktop/XMAS_analysis/outputs/s06_MS_XMAS_relabeled.rds')


# In[7]:


cellinfo <- sec@meta.data[,c('seurat_clusters_BiMod_OC','bimod_oc','bimod_oc_des','orig.ident')]
colnames(cellinfo) <- c('seurat_clusters_BiMod_OC','bimod_oc','bimod_oc_des','orig.ident')


# ## 4.1 RSS in bimod_oc_des 

# In[8]:


cellTypes <-  as.data.frame(subset(cellinfo,select = 'bimod_oc_des'))
selectedResolution <- "bimod_oc_des"
cellAnnotation = cellTypes[colnames(regulonAUC),
                           selectedResolution]
cellAnnotation = na.omit(cellAnnotation)


# In[9]:


rss <- calcRSS(AUC = getAUC(regulonAUC),
               cellAnnotation = cellAnnotation)
rss = na.omit(rss)


# In[10]:


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


# In[11]:


rss


# In[17]:


Idents(sec) <- 'bimod_oc_des'


# In[20]:


options(repr.plot.width=10, repr.plot.height=20)

plot_list <- list()

for (i in seq(1, 15)) {
  setName <- table(sec@active.ident)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 5) 


# In[21]:


rss <- rss[!rownames(rss) %in% c("BRCA1(+)","MYC(+)","HOXB5(+)","THRB(+)"),]


# In[23]:


options(repr.plot.width=10, repr.plot.height=20)

plot_list <- list()

for (i in seq(1, 15)) {
  setName <- table(sec$bimod_oc_des)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 5) 


# In[24]:


options(repr.plot.width=3.75, repr.plot.height=4.5)

plot_list <- list()
for (i in seq(1, 15)) {
  setName <- table(sec$bimod_oc_des)[i]%>%names()
  print(plotRSS_oneSet(rss[, ], setName = setName, n = 5))
}


# ## 4.2 RSS in bimod_oc

# In[44]:


cellTypes <-  as.data.frame(subset(cellinfo,select = 'bimod_oc'))
selectedResolution <- "bimod_oc"
cellAnnotation = cellTypes[colnames(regulonAUC),
                           selectedResolution]
cellAnnotation = na.omit(cellAnnotation)


# In[45]:


rss <- calcRSS(AUC = getAUC(regulonAUC),
               cellAnnotation = cellAnnotation)
rss = na.omit(rss)


# In[46]:


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


# In[47]:


rss


# In[48]:


Idents(sec) <- 'bimod_oc'


# In[51]:


options(repr.plot.width=10, repr.plot.height=10)

plot_list <- list()

for (i in seq(1, 8)) {
  setName <- table(sec@active.ident)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 2) 


# In[52]:


rss <- rss[!rownames(rss) %in% c("BRCA1(+)","MYC(+)","HOXB5(+)","THRB(+)"),]


# In[58]:


options(repr.plot.width=10, repr.plot.height=10)

plot_list <- list()

for (i in seq(1, 8)) {
  setName <- table(sec$bimod_oc)[i]%>%names()
  plot_list[[i]] <- plotRSS_oneSet(rss[, ], setName = setName, n = 5)
}

# Arrange plots in a grid
plot_grid(plotlist = plot_list, nrow = 2) 


# In[59]:


options(repr.plot.width=3.75, repr.plot.height=4.5)

plot_list <- list()
for (i in seq(1, 8)) {
  setName <- table(sec$bimod_oc)[i]%>%names()
  print(plotRSS_oneSet(rss[, ], setName = setName, n = 5))
}


# ## 4.3 Regulon heatmap 

# In[64]:


library(ggheatmap)
library(reshape2)
library(RColorBrewer)

#regulon表达水平热图
tfs <- c("Nr2f1(+)","Cebpd(+)","Hnf4g(+)","Twist1(+)","Twist2(+)","Prrx2(+)",
         "Lef1(+)","Foxl2(+)","Foxp1(+)","Hey2(+)","Sox6(+)","Msx1(+)")
rss_data <- rssPlot$plot$data[which(rssPlot$plot$data$Topic %in% tfs),]
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
col_ann <- data.frame(group= c(rep("Acta2+ SMC",1),
                               rep("Notch3+ FB",1),
                               rep("Col8a1+ FB",1),
                               rep("Kcnma1+ SMC",1),
                               rep("Kdr+ EC",1),
                               rep("Ednrb+ EC",1),
                               rep("Pecam1+ EC",1),
                               rep("Pi16+ FB",1),
                               rep("Eln+ FB",1),
                               rep("Ly6c1+ EC",1),
                               rep("Pdgfra+ FB",1),
                               rep("Klf4+ EC",1),
                               rep("Ednra+ SMC",1),
                               rep("Angpt1+ FB",1)))
rownames(col_ann) <- colnames(rss_data)
groupcol <- colorRampPalette(brewer.pal(14,'Set3'))(14)
names(groupcol) <- c("Acta2+ SMC","Notch3+ FB","Col8a1+ FB","Kcnma1+ SMC","Kdr+ EC",    
                     "Ednrb+ EC","Pecam1+ EC","Pi16+ FB","Eln+ FB","Ly6c1+ EC",  
                     "Pdgfra+ FB","Klf4+ EC","Ednra+ SMC","Angpt1+ FB")
col <- list(group=groupcol)
text_columns <- sample(colnames(rss_data),0)
ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns) #转录因子表达水平热图
top3tfgene <- c("Nr2f1","Cebpd","Hnf4g","Twist1","Twist2","Prrx2",
                "Lef1","Foxl2","Foxp1","Hey2","Sox6","Msx1")
top3gene_cell_exp <- AverageExpression(sc,
                                       assays = 'RNA',
                                       features = top3tfgene,
                                       group.by = 'celltype',
                                       slot = 'data') 
top3gene_cell_exp <- as.data.frame(top3gene_cell_exp$RNA)
top3marker_exp <- t(scale(t(top3gene_cell_exp),scale = T,center = T))
ggheatmap(top3marker_exp,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
          cluster_rows = T,cluster_cols = F,scale = "row",
          annotation_cols = col_ann,
          annotation_color = col,
          legendName="Relative value",
          text_show_cols = text_columns)
# In[ ]:





# In[ ]:




