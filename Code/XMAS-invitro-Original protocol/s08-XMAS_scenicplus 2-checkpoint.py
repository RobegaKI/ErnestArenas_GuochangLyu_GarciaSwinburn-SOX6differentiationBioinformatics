#!/usr/bin/env python
# coding: utf-8

# # Preprocessing

# In[1]:


library(Seurat)
library(SeuratDisk)


# In[3]:


sec <- readRDS('~/Desktop/XMAS_analysis/outputs/s06_MS_XMAS_relabeled.rds')


# In[6]:


sec


# In[7]:


sec[['peaks']] <- NULL
sec[['GenePromAcc']] <- NULL
sec[['unspliced']] <- NULL


# In[8]:


sec


# In[11]:


head(sec@assays$RNA@scale.data)


# In[9]:


head(sec@assays$RNA@data)


# In[10]:


head(sec@assays$RNA@counts)


# In[12]:


SaveH5Seurat(sec, filename = ("~/Desktop/XMAS_analysis/outputs/s08_XMAS_scplus_raw_included_no_updated_metadata.h5Seurat"))
SeuratDisk::Convert("~/Desktop/XMAS_analysis/outputs/s08_XMAS_scplus_raw_included_no_updated_metadata.h5Seurat", dest = "h5ad")


# In[1]:


import anndata
import scanpy as sc


# In[2]:


adata = sc.read_h5ad('/Users/gclyu07/Desktop/XMAS_analysis/outputs/s08_XMAS_scplus_raw_included_no_updated_metadata.h5ad')


# In[3]:


adata.obs['celltype.split'] = adata.obs['bimod_oc'].astype(str) + " " + adata.obs['orig.ident'].astype(str)
adata.obs['celltype.split'] = adata.obs['celltype.split'].astype('category')


# In[4]:


adata.obs['celltype.split.2'] = adata.obs['bimod_oc_des'].astype(str) + " " + adata.obs['orig.ident'].astype(str)
adata.obs['celltype.split.2'] = adata.obs['celltype.split.2'].astype('category')


# In[7]:


adata


# In[9]:


adata.raw = adata


# In[10]:


adata.raw.X


# In[11]:


import pycisTopic
import pickle


# In[12]:


with open('/Users/gclyu07/Desktop/XMAS_analysis/outputs/cistopic_obj.pkl', "rb") as file:
    cistopic_obj = pickle.load(file)


# In[13]:


adata_cis = adata[adata.obs_names.isin(cistopic_obj.cell_names)]


# In[15]:


len(cistopic_obj.cell_names)


# In[14]:


adata_cis


# In[17]:


adata_cis.write_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/s08_adata_pycisTopic_matched.h5ad")


# # Scenicplus 

# In[1]:


import os
os.chdir("/Users/gclyu07/Desktop/XMAS_analysis/outputs/scplus_outputs/")


# In[93]:


import mudata
import scenicplus


# In[182]:


scplus_mdata = mudata.read("./scplusmdata.h5mu")


# In[4]:


scplus_mdata.uns["direct_e_regulon_metadata"]


# In[5]:


scplus_mdata.uns["extended_e_regulon_metadata"]


# In[58]:


eRegulon_gene_AUC.obs_names


# In[77]:


from scenicplus.scenicplus_class import mudata_to_scenicplus


# In[200]:


scplus_obj = mudata_to_scenicplus(
    mdata = scplus_mdata,
    path_to_cistarget_h5 = "./ctx_results.hdf5",
    path_to_dem_h5 = "./dem_results.hdf5"
)


# In[82]:


scplus_obj.uns


# In[110]:


import pickle
pickle.dump(
    scplus_obj,
    open("/Users/gclyu07/Desktop/XMAS_analysis/outputs/scplus_outputs/scplus_obj.pkl", "wb"))


# ## eRegulon dimensionality reduction

# In[183]:


import scanpy as sc
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)


# In[184]:


eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]


# In[185]:


sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")


# In[186]:


sc.tl.umap(eRegulon_gene_AUC)


# In[10]:


sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:bimod_oc_des")


# In[12]:


sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:orig.ident")


# In[25]:


sc.set_figure_params(figsize=(7, 5))
sc.pl.umap(eRegulon_gene_AUC, color=['scRNA_counts:bimod_oc_des'], add_outline=True, size=25,
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Cell types'])


# In[26]:


sc.set_figure_params(figsize=(7, 5))
oi_colors={ "D11": "#F0AD4E", "D16": "#D9534F", "D28": "#428BCA", "D42": "#9933CC", "D56":"#66CCCC" }
sc.pl.umap(eRegulon_gene_AUC, color=['scRNA_counts:orig.ident'], add_outline=True, size=25, palette=oi_colors,
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Identities'])


# ## eRegulon specificity score

# In[187]:


from scenicplus.RSS import (regulon_specificity_scores, plot_rss)


# In[188]:


rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:bimod_oc_des",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)


# In[47]:


rss


# In[54]:


plot_rss(
    data_matrix = rss,
    top_n = 5,
    num_columns = 5,
    figsize = (7, 8),
)


# ## Plot eRegulon enrichment scores

# In[31]:


sc.pl.umap(eRegulon_gene_AUC, color = list(set([x for xs in [rss.loc[ct].sort_values()[0:2].index for ct in rss.index] for x in xs ])))


# In[62]:


rss


# In[69]:


sc.set_figure_params(figsize=(5, 4))
sc.pl.umap(eRegulon_gene_AUC, color = 'STAT4_direct_+/+_(10g)', cmap='viridis_r')


# ## Heatmap dotplot

# In[41]:


from scenicplus.plotting.dotplot import heatmap_dotplot


# In[46]:


heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:bimod_oc_des",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    figsize = (70, 8)
)


# # Downstream analysis, export and plotting

# In[202]:


scenicplus.eregulon_enrichment.binarize_AUC(scplus_obj,
                                           auc_key = 'eRegulon_AUC',
                                           out_key = 'eRegulon_AUC_thresholds',
                                           signature_keys = ['Gene_based', 'Region_based'],
                                           n_cpu = 1)


# In[220]:


# Format eRegulons
from scenicplus.eregulon_enrichment import *
get_eRegulons_as_signatures(scplus_obj.uns['eRegulon_metadata'])


# In[230]:


scplus_obj.uns


# In[ ]:


scenicplus.eregulon_enrichment.rank_data(scplus_obj.uns['eRegulon_metadata'], axis = 1, seed = 1121)

