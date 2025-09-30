#!/usr/bin/env python
# coding: utf-8

# # 1. Celltypist

# In[1]:


import scanpy as sc
import scvelo as scv
import celltypist
import time
import numpy as np


# In[2]:


mid = sc.read_h5ad('/Users/gclyu07/Desktop/R_data/Kamath-CellTypist.h5ad',)


# In[206]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(10, 5))
sc.pl.umap(mid, color = ['Cell_Type'], legend_loc = 'on data',add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25,)


# In[208]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(mid, color=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'],
           title=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'], ncols=4,
           cmap='viridis_r')


# ## 1.1 DA

# In[69]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/R_data/outputs/0326_DA_annotated.h5ad")


# In[70]:


adata


# In[71]:


adata.obs['seurat_clusters_pca_RNA'] = adata.obs['seurat_clusters_pca_RNA'].astype(str)
adata.obs['seurat_clusters_pca_RNA_2'] = adata.obs['seurat_clusters_pca_RNA_2'].astype(str)
adata.obs['seurat_clusters_pca_RNA_3'] = adata.obs['seurat_clusters_pca_RNA_3'].astype(str)


# In[72]:


adata.obs['orig.ident'].value_counts()


# In[7]:


adata.obsm['X_umap'] = adata.obsm['X_umap']


# In[175]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2'],frameon=False)


# In[50]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'],
           title=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'], ncols=4, use_raw=True,
           cmap='viridis_r')


# In[14]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["TH","SOX6"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='TH_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='TH_score', cmap='viridis_r')
sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["TH","CALB1"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='TH_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='TH_score', cmap='viridis_r')


# In[61]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["TH","SOX6"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='TH_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='TH_score', cmap='viridis_r')


# In[15]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["OTX2", "LMX1A", "EN1", "SIM2","FOXA2"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='mid_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='mid_score', cmap='viridis_r')
sc.pl.violin(
        adata,
        ["mid_score"],
        groupby="orig.ident",
        stripplot=False,  # remove the internal dots
        inner="box",  # adds a boxplot inside violins
    )


# In[176]:


sc.pl.stacked_violin(
    adata, ["TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"],
    groupby="seurat_clusters_pca_RNA_2", swap_axes=False, dendrogram=True
)


# In[ ]:





# In[29]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 1.5))
sc.pl.umap(adata, color = ['orig.ident','seurat_clusters_pca_RNA_2'], legend_loc = 'on data',add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25,)


# In[ ]:





# In[20]:


sc.pp.normalize_total(mid, target_sum = 1e4)
sc.pp.log1p(mid)


# In[73]:


adata.X


# In[30]:


mid


# In[141]:


# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
t_start = time.time()
model_fs = celltypist.train(mid, 'Cell_Type', n_jobs = 10, max_iter = 20, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[142]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -500, axis = 1)[:, -500:]


# In[143]:


gene_index = np.unique(gene_index)


# In[144]:


print(f"Number of genes selected: {len(gene_index)}")


# In[145]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid, 'Cell_Type', check_expression = False, n_jobs = 10, max_iter = 400)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[146]:


# Save the model.
model.write('/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model.pkl')


# In[147]:


adata


# In[148]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[149]:


adata = predictions.to_adata()


# In[150]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2', 'predicted_labels','majority_voting'])


# In[151]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'majority_voting')


# In[ ]:





# In[117]:


import numpy as np

mid.obs['group'] = np.where(
    mid.obs['Cell_Type'].str.startswith('CALB1'), 
    'CALB1', 
    np.where(
        mid.obs['Cell_Type'].str.startswith('SOX6'), 
        'SOX6', 
        'Other'  # Optional: for non-matching clusters
    )
)

print(mid.obs['group'].value_counts())


# In[152]:


# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
t_start = time.time()
model_fs = celltypist.train(mid, 'group', n_jobs = 10, max_iter = 20, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[162]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -1000, axis = 1)[:, -1000:]


# In[163]:


gene_index = np.unique(gene_index)


# In[164]:


print(f"Number of genes selected: {len(gene_index)}")


# In[166]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid, 'group', check_expression = False, n_jobs = 10, max_iter = 400)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[167]:


# Save the model.
model.write('/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_2.pkl')


# In[168]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_2.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[169]:


adata = predictions.to_adata()


# In[170]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2', 'predicted_labels','majority_voting'])


# In[171]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'majority_voting')


# In[ ]:





# ## 1.2 Neu

# In[177]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/R_data/outputs/0325_neu_annotated.h5ad")


# In[178]:


adata


# In[179]:


adata.obs['seurat_clusters_pca_RNA'] = adata.obs['seurat_clusters_pca_RNA'].astype(str)
adata.obs['seurat_clusters_pca_RNA_2'] = adata.obs['seurat_clusters_pca_RNA_2'].astype(str)


# In[180]:


adata.obs['orig.ident'].value_counts()


# In[181]:


adata.obsm['X_umap'] = adata.obsm['X_umap']


# In[182]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2'],frameon=False)


# In[183]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'],
           title=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'], ncols=4, use_raw=True,
           cmap='viridis_r')


# In[184]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["TH","SOX6"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='TH_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='TH_score', cmap='viridis_r')
sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["TH","CALB1"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='TH_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='TH_score', cmap='viridis_r')


# In[186]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["OTX2", "LMX1A", "EN1", "SIM2","FOXA2"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='mid_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='mid_score', cmap='viridis_r')
sc.pl.violin(
        adata,
        ["mid_score"],
        groupby="orig.ident",
        stripplot=False,  # remove the internal dots
        inner="box",  # adds a boxplot inside violins
    )


# In[189]:


adata


# In[187]:


sc.pl.stacked_violin(
    adata, ["TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"],
    groupby="seurat_clusters_pca_RNA_2", swap_axes=False, dendrogram=True
)


# In[ ]:





# In[195]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(5, 3))
sc.pl.umap(adata, color = ['orig.ident','seurat_clusters_pca_RNA_2'], legend_loc = 'on data',add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=2,frameon=False,  size=25,)


# In[197]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[198]:


adata = predictions.to_adata()


# In[199]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2', 'predicted_labels','majority_voting'])


# In[200]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'majority_voting')


# In[201]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_2.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[202]:


adata = predictions.to_adata()


# In[203]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2', 'predicted_labels','majority_voting'])


# In[204]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'majority_voting')


# In[ ]:





# ## 1.3 DA-refined

# In[229]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/R_data/outputs/0326_DA_annotated.h5ad")


# In[230]:


adata


# In[231]:


adata.obs['seurat_clusters_pca_RNA'] = adata.obs['seurat_clusters_pca_RNA'].astype(str)
adata.obs['seurat_clusters_pca_RNA_2'] = adata.obs['seurat_clusters_pca_RNA_2'].astype(str)
adata.obs['seurat_clusters_pca_RNA_3'] = adata.obs['seurat_clusters_pca_RNA_3'].astype(str)


# In[232]:


adata.obs['orig.ident'].value_counts()


# In[233]:


adata.obsm['X_umap'] = adata.obsm['X_umap']


# In[253]:


# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
t_start = time.time()
model_fs = celltypist.train(mid, 'Cell_Type', n_jobs = 10, max_iter = 20)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[254]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -500, axis = 1)[:, -500:]


# In[255]:


gene_index = np.unique(gene_index)


# In[256]:


print(f"Number of genes selected: {len(gene_index)}")


# In[257]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid, 'Cell_Type', n_jobs = 10, max_iter = 400)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[258]:


# Save the model.
model.write('/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_3.pkl')


# In[259]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_3.pkl',
                                  majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[260]:


adata = predictions.to_adata()


# In[261]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2', 'predicted_labels','majority_voting'])


# In[262]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA_3', use_as_prediction = 'predicted_labels')


# In[ ]:





# In[117]:


import numpy as np

mid.obs['group'] = np.where(
    mid.obs['Cell_Type'].str.startswith('CALB1'), 
    'CALB1', 
    np.where(
        mid.obs['Cell_Type'].str.startswith('SOX6'), 
        'SOX6', 
        'Other'  # Optional: for non-matching clusters
    )
)

print(mid.obs['group'].value_counts())


# In[152]:


# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
t_start = time.time()
model_fs = celltypist.train(mid, 'group', n_jobs = 10, max_iter = 20, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[162]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -1000, axis = 1)[:, -1000:]


# In[163]:


gene_index = np.unique(gene_index)


# In[164]:


print(f"Number of genes selected: {len(gene_index)}")


# In[166]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid, 'group', check_expression = False, n_jobs = 10, max_iter = 400)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[167]:


# Save the model.
model.write('/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_2.pkl')


# In[168]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/R_data/Kamath-CellTypist_model_2.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[169]:


adata = predictions.to_adata()


# In[170]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2', 'predicted_labels','majority_voting'])


# In[171]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA', use_as_prediction = 'majority_voting')


# # 2. MetaNeighbor

# In[152]:


mid


# In[6]:


library(MetaNeighbor)
library(SummarizedExperiment)


# In[292]:


DA@meta.data$Dataset <- 'Graft'
mid$Dataset <- "Reference"


# In[163]:


Idents(mid)=mid@meta.data$Cell_Type
Idents(DA)=DA@meta.data$seurat_clusters_pca_RNA_2


# In[164]:


ob.list <- list(mid,DA)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)
mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)


# In[294]:


DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)

table(Idents(object = mda.integrated),mda.integrated@meta.data$Dataset)
var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[166]:


Study_ID = rep(c('1', '2'), c(ncol(mid),ncol(DA)))
Celltype = c(as.character(mid@meta.data$Cell_Type),as.character(DA$seurat_clusters_pca_RNA_2))


# In[167]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[168]:


length(unique(mid@meta.data$Cell_Type))
length(unique(DA$seurat_clusters_pca_RNA_2))


# In[190]:


table(DA$seurat_clusters_pca_RNA_2)


# In[174]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 10), rep('Experiment',4)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=8, repr.plot.height=8)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "average",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[175]:


Study_ID = rep(c('1', '2'), c(ncol(mid),ncol(DA)))
Celltype = c(as.character(mid@meta.data$Cell_Type),as.character(DA$seurat_clusters_pca_RNA_3))


# In[176]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[177]:


length(unique(mid@meta.data$Cell_Type))
length(unique(DA$seurat_clusters_pca_RNA_3))


# In[178]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 10), rep('Experiment',12)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=8, repr.plot.height=8)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "average",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[183]:


Study_ID = rep(c('1', '2'), c(ncol(mid),ncol(DA)))
Celltype = c(as.character(mid@meta.data$Cell_Type),as.character(DA$seurat_clusters_pca_RNA))


# In[184]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[185]:


length(unique(mid@meta.data$Cell_Type))
length(unique(DA$seurat_clusters_pca_RNA))


# In[189]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 10), rep('Experiment',4)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=8, repr.plot.height=8)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "median",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[ ]:





# # 3. Dotplots

# ## da marker genes

# In[192]:


options(repr.plot.width=10, repr.plot.height=4)
plot <- DotPlot(DA, assay = "RNA", dot.min = 0.03, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "NEUROG2", "DCX"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[195]:


options(repr.plot.width=10, repr.plot.height=4)
plot <- DotPlot(DA, assay = "RNA", group.by = 'seurat_clusters_pca_RNA', dot.min = 0.01, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6","KCNJ6","CALB1"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# ## cluster marker genes

# In[199]:


table(DA$seurat_clusters_pca_RNA_2)


# In[266]:


DefaultAssay(DA) <- "RNA"
Idents(object = DA) <- "seurat_clusters_pca_RNA_2"

markers <- FindAllMarkers(DA, slot = "data", min.pct = 0.2, logfc.threshold = 0.5, only.pos = TRUE)
markers_best <- markers[markers$p_val_adj < 0.05,]


# In[267]:


markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7
top7


# In[ ]:





# In[202]:


DefaultAssay(DA) <- "RNA"
Idents(object = DA) <- "seurat_clusters_pca_RNA_2"

markers <- FindMarkers(DA, slot = "data", min.pct = 0.1, ident.1 = c(7,10), ident.2 = c(1,13), logfc.threshold = 0.5, only.pos = TRUE)
markers_best <- markers[markers$p_val_adj < 0.05,]


# In[204]:


markers_best


# In[210]:


markers_best %>% top_n(n = 20, wt = avg_log2FC) -> top20


# In[217]:


top20


# In[225]:


gene_names <- rownames(top20)
filtered_genes <- gene_names[!grepl("^ENSG", gene_names)]


# In[254]:


filtered_genes


# In[231]:


options(repr.plot.width=15, repr.plot.height=4)
plot <- DotPlot(DA, assay = "RNA", group.by = 'seurat_clusters_pca_RNA_2', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = unique(filtered_genes))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[ ]:





# In[31]:


# Create a new metadata column for merged clusters
DA$merged_clusters <- as.character(DA$seurat_clusters_pca_RNA_2)
# Merge clusters 7 and 10 into a new group (e.g., "7_10")
DA$merged_clusters[DA$seurat_clusters_pca_RNA_2 %in% c("7", "10")] <- "7_10"


# In[32]:


DA <- subset(DA, subset = seurat_clusters_pca_RNA_2 != "0")


# In[33]:


table(DA$merged_clusters)


# In[278]:


options(repr.plot.width=15, repr.plot.height=4)
VlnPlot(DA, group.by = 'merged_clusters', features = c("SOX6", "PITX3", "ONECUT2", "KCNJ6","CALB1", "VIP"), ncol = 6)


# In[35]:


DA$merged_clusters <- factor(DA$merged_clusters)
levels(DA$merged_clusters) <- c("DA_CALB1", "DA_ONECUT2", "DA_SOX6")


# In[36]:


table(DA$merged_clusters)


# ______

# # 4. Volcano plot

# In[232]:


library(scRNAtoolVis)


# In[233]:


high_contrast_colors <- colorspace::rainbow_hcl(11)


# In[234]:


Idents(object = DA) <- "seurat_clusters_pca_RNA_2"
deg <- FindAllMarkers(DA, min.pct = 0.1, logfc.threshold = 0.5, only.pos = FALSE)


# In[241]:


options(repr.plot.width=10, repr.plot.height=5)
Idents(AMS) <- "seurat_clusters_pca_RNA_2"
jjVolcano(diffData = deg %>% filter(!grepl("^RP|^ENSG", gene)),
          tile.col = high_contrast_colors,
          pvalue.cutoff = 0.01) + NoLegend()


# In[248]:


options(repr.plot.width=8, repr.plot.height=7)
markerVocalno(markers = deg %>% filter(!grepl("^RP|^ENSG", gene)),labelCol = ggsci::pal_npg()(11)) + NoLegend()


# In[ ]:





# In[ ]:





# In[281]:


Idents(object = DA) <- "merged_clusters"
deg <- FindAllMarkers(DA, min.pct = 0.1, logfc.threshold = 0.5, only.pos = FALSE)


# In[283]:


options(repr.plot.width=6, repr.plot.height=5)
Idents(AMS) <- "merged_clusters"
jjVolcano(diffData = deg %>% filter(!grepl("^RP|^ENSG", gene)),
          tile.col = high_contrast_colors,
          pvalue.cutoff = 0.01) + NoLegend()


# In[285]:


options(repr.plot.width=6, repr.plot.height=7)
markerVocalno(markers = deg %>% filter(!grepl("^RP|^ENSG", gene)),labelCol = ggsci::pal_npg()(11)) + NoLegend()


# In[ ]:





# ### METANEIGHBOR

# In[328]:


Study_ID = rep(c('1', '2'), c(ncol(mid),ncol(DA)))
Celltype = c(as.character(mid@meta.data$Cell_Type),as.character(DA$merged_clusters))


# In[329]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[330]:


length(unique(mid@meta.data$Cell_Type))
length(unique(DA$merged_clusters))


# In[331]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(Dataset=c(rep("Reference", 10), rep('Experiment',3)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=10, repr.plot.height=10)
pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "average",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[362]:


Study_ID = rep(c(' ', '  '), c(ncol(mid),ncol(DA)))
Celltype = c(as.character(mid@meta.data$Cell_Type),as.character(DA$merged_clusters))


# In[363]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[364]:


rownames(celltype_NV) <- gsub("\\|", "", rownames(celltype_NV))
colnames(celltype_NV) <- gsub("\\|", "", colnames(celltype_NV))


# In[365]:


length(unique(mid@meta.data$Cell_Type))
length(unique(DA$merged_clusters))


# In[366]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(Dataset=c(rep("Reference", 10), rep('Experiment',3)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=12, repr.plot.height=10)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "average",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[225]:


options(repr.plot.width=10, repr.plot.height=4)
VlnPlot(DA, group.by = 'merged_clusters', features = c("SOX6", "ONECUT2", "CALB1","NR4A2"), ncol = 4)


# In[375]:


options(repr.plot.width=6, repr.plot.height=8)
VlnPlot(DA, group.by = 'merged_clusters', features = c("SOX6", "ONECUT2", "CALB1","NR4A2"), ncol = 2)


# In[236]:


modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


# In[245]:


options(repr.plot.width=15, repr.plot.height=10)
StackedVlnPlot(mid, group.by = 'Cell_Type', features = c("CNTNAP4","SLIT2","DSCAML1","NFIA","SOX6","CALB1"))


# In[ ]:





# In[40]:


devtools::install_github('kevinblighe/EnhancedVolcano')


# In[41]:


library(EnhancedVolcano)


# In[43]:


Idents(DA) <- 'merged_clusters'
deg <- FindMarkers(DA, slot = "data", min.pct = 0.1,logfc.threshold = 0.5,
                   ident.1 = c("DA_SOX6"), ident.2 = c("DA_ONECUT2", "DA_CALB1"))


# In[63]:


options(repr.plot.width=5, repr.plot.height=4)
VlnPlot(DA, group.by = 'merged_clusters', features = c("CNTNAP4"))


# In[247]:


options(repr.plot.width=5, repr.plot.height=4)
VlnPlot(DA, group.by = 'merged_clusters', features = c("ZNF536"))


# In[49]:


EnhancedVolcano(deg,
    lab = rownames(deg),
    x = 'avg_log2FC',
    y = 'p_val')


# In[60]:


options(repr.plot.width=6, repr.plot.height=8)

EnhancedVolcano(deg,
    lab = rownames(deg),
    x = 'avg_log2FC',
    y = 'p_val',
    pCutoff = 10e-32,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 6.0)


# In[96]:


deg


# In[145]:


deg <- deg %>% filter(!grepl("^RP|^ENSG", rownames(deg)))
deg$gene <- rownames(deg)
markers_best <- deg[deg$p_val_adj < 0.05,]


# In[148]:


markers_best[1:10,]

selected_balanced <- markers_best %>%
  group_by(sign(avg_log2FC)) %>%      # Split by positive/negative FC
  arrange(desc(abs(avg_log2FC))) %>%  # Sort within each group
  slice_head(n = 5) %>%               # Take top 5 from each
  ungroup()
# In[160]:


markers_best[1:10,]$gene


# In[246]:


options(repr.plot.width=7, repr.plot.height=8)
  EnhancedVolcano(deg,
    lab = rownames(deg),
    x = 'avg_log2FC',
    y = 'p_val',
    #selectLab = c(markers_best[1:10,]$gene,"SOX6"),
    selectLab = c("CNTNAP4","NFIA","SLIT2","DSCAML1","SOX6","KLHL1","MGAT4C","ZNF536"),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 1.0,
    pointSize = 3.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    parseLabels = TRUE,
    col = c('darkgrey', 'darkgrey', 'blue', 'red3'),
    colAlpha = 4/5,
    legendPosition = 'bottom',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black') + coord_flip()


# In[180]:


table(DA@meta.data$orig.ident,DA@meta.data$merged_clusters)


# In[204]:


desired_order <- c("DA_SOX6", "DA_CALB1", "DA_ONECUT2")


# In[205]:


table_samples_by_cell_type <- DA@meta.data %>%
  dplyr::group_by(orig.ident, merged_clusters) %>%
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(merged_clusters, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(dplyr::select(., all_of(desired_order)))) %>%
  dplyr::select(c("orig.ident", "total_cell_count", all_of(desired_order)))


# In[224]:


table_samples_by_cell_type


# In[214]:


options(repr.plot.width=4.5, repr.plot.height=7)

temp_labels <- DA@meta.data %>%
  group_by(orig.ident) %>%
  tally()


table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(desired_order))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(
      x = orig.ident,
      y = Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ),
    color = 'black', size = 2.8
  ) +

  scale_fill_manual(name = 'Neuron type', values = c('#E63863', '#DCC1DD', '#CCE0F5')) +
  scale_y_continuous(
    name = 'Percentage [%]',
    labels = scales::percent_format(),
    expand = c(0.01,0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  ) 


# In[223]:


options(repr.plot.width=8, repr.plot.height=2)

temp_labels <- DA@meta.data %>%
  group_by(orig.ident) %>%
  tally()


table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(desired_order))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(
      x = orig.ident,
      y = Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ),
    color = 'black', size = 0
  ) +

  scale_fill_manual(name = 'Neuron type', values = c('#E63863', '#DCC1DD', '#CCE0F5')) +
  scale_y_continuous(
    name = 'Percentage [%]',
    labels = scales::percent_format(),
    expand = c(0.01,0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  )  + coord_flip()


# # 5. ClusterProfiler

# In[376]:


Idents(DA) <- 'merged_clusters'
for (ct in unique(DA$merged_clusters)) {
    
cluster.markers <- FindMarkers(DA, ident.1 = ct,min.pct = 0.25)  

#RP,ENSG gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP|^ENSG", rownames(cluster.markers)))
    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  
options(repr.plot.width=7.2, repr.plot.height=3)
print(ego.bp %>% filter(p.adjust < 0.03) %>%
ggplot(showCategory = 3,
  aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab(NULL) +
  ylab(NULL) + 
  ggtitle(paste0(ct)) +
  scale_y_discrete(labels = label_wrap(45)) + NoLegend() +
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}


# In[377]:


Idents(DA) <- 'merged_clusters'
for (ct in unique(DA$merged_clusters)) {
    
cluster.markers <- FindMarkers(DA, ident.1 = ct, min.pct = 0.5, only.pos = TRUE)  

#RP,ENSG gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP|^ENSG", rownames(cluster.markers)))
    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  
options(repr.plot.width=7.2, repr.plot.height=3)
print(ego.bp %>% filter(p.adjust < 0.03) %>%
ggplot(showCategory = 3,
  aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab(NULL) +
  ylab(NULL) + 
  ggtitle(paste0(ct)) +
  scale_y_discrete(labels = label_wrap(45)) + NoLegend() +
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[395]:


Idents(DA) <- 'merged_clusters'
for (ct in unique(DA$merged_clusters)) {
    
cluster.markers <- FindMarkers(DA, ident.1 = ct,min.pct = 0.1, logfc.threshold = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP|^ENSG", rownames(cluster.markers)))
    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  
options(repr.plot.width=7, repr.plot.height=3.5)        
print(ego.bp %>% filter(p.adjust < 0.03) %>%
ggplot(showCategory = 5,
  aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab(NULL) +
  ylab(NULL) + 
  ggtitle(paste0(ct)) +
  scale_y_discrete(labels = label_wrap(30)) + 
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}


# # 7. Python illustration

# In[211]:


library(Seurat)
library(SeuratDisk)


# In[396]:


SaveH5Seurat(DA, filename = paste0(OS_path_outputs,"0401_DA_Merged.h5Seurat"))


# In[397]:


SeuratDisk::Convert(paste0(OS_path_outputs,"0401_DA_Merged.h5Seurat"), dest = "h5ad")


# In[2]:


import scanpy as sc
import scvelo as scv
import celltypist
import time
import numpy as np


# In[3]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/R_data/outputs/0401_DA_Merged.h5ad")


# In[4]:


adata


# In[5]:


adata.obs['merged_clusters'].value_counts()


# In[6]:


alist = {
    0: "DA_CALB1", 
    1: "DA_ONECUT2", 
    2: "DA_SOX6",
}
adata.obs['merged_clusters'] = (
adata.obs['merged_clusters']
.astype('category')
.map(alist)
)


# In[7]:


adata.obs['merged_clusters'].value_counts()


# In[8]:


adata.obsm['X_umap'] = adata.obsm['X_umap']


# In[8]:


sc.pl.umap(adata, color = ['merged_clusters'])


# In[24]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 2))
sc.pl.umap(adata, color = ['orig.ident','merged_clusters'],title=['Original identity','DA subtypes'],
           add_outline=True,  legend_fontsize=6, legend_fontoutline=1,frameon=False, size=50,)


# In[13]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 2))
sc.pl.umap(adata, color = ['merged_clusters'],title=['DA subtypes'],
           add_outline=True,  legend_fontsize=6, legend_fontoutline=1,frameon=False, size=50,)


# In[22]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 2))
sc.pl.umap(adata, color=['CALCRL','ADAMTS6'], ncols=4, use_raw=True,
           cmap='viridis_r')


# In[35]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(3, 2))
sc.pl.umap(adata, color=['MMP2'], ncols=4, use_raw=True,
           cmap='viridis_r')


# In[ ]:





# In[ ]:





# In[36]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['LMX1A', 'EN1','OTX2','TH','DDC','NR4A2','SLC18A2','SLC6A3','SOX6','CALB1','ONECUT2','PITX3'],
           title=['LMX1A', 'EN1','OTX2','TH','DDC','NR4A2','SLC18A2','SLC6A3','SOX6','CALB1','ONECUT2','PITX3'], ncols=4, use_raw=True,
           cmap='viridis_r')


# In[31]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['LMX1A', 'EN1','OTX2','TH','DDC','NR4A2','SLC18A2','SLC6A3','PITX3','SOX6','CALB1','ONECUT2'],
           title=['LMX1A', 'EN1','OTX2','TH','DDC','NR4A2','SLC18A2','SLC6A3','PITX3','SOX6','CALB1','ONECUT2'], ncols=3, use_raw=True,
           cmap='viridis_r')


# In[37]:


adata.write("/Users/gclyu07/Desktop/R_data/outputs/0409_DA.h5ad")


# In[40]:


adata.write_loom("/Users/gclyu07/Desktop/R_data/outputs/0409_DA.loom", write_obsm_varm=True)

