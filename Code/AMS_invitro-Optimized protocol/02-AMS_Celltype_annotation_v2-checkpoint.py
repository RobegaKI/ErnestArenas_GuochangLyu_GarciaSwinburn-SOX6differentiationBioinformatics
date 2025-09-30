#!/usr/bin/env python
# coding: utf-8

# # 1. Celltypist

# In[1]:


library(Seurat)
library(SeuratDisk)


# In[201]:


SaveH5Seurat(AMS, filename = paste0(OS_path_outputs,"0725_AMS_Celltypist.h5Seurat"))


# In[202]:


SeuratDisk::Convert(paste0(OS_path_outputs,"0725_AMS_Celltypist.h5Seurat"), dest = "h5ad")


# In[ ]:





# In[1]:


import scanpy as sc
import scvelo as scv
import celltypist
import time
import numpy as np

mid = sc.read_h5ad('Desktop/XMAS_analysis/outputs_no_d0/SL_midbrain_2.h5ad',)mid_2 = sc.read_h5ad('Desktop/XMAS_analysis/outputs/f01-XMAS_oc_monocle_unspliced.h5ad')
# In[2]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/AMS_analysis/outputs/0725_AMS_Celltypist.h5ad")


# In[4]:


adata.obs['seurat_clusters_pca_RNA_0.3'] = adata.obs['seurat_clusters_pca_RNA_0.3'].astype(str)


# In[5]:


adata.obs['orig.ident'].value_counts()


# In[20]:


adata.obsm['X_umap'] = adata.obsm['X_umap_pca_RNA_0.3']


# In[21]:


sc.pl.umap(adata, color = ['Cell_type_SL_abv0.3'],frameon=False)


# In[22]:


sc.pl.umap(adata, color = ['orig.ident','seurat_clusters_pca_RNA_0.3'], frameon=False)


# In[23]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'],
           title=['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'], ncols=4, use_raw=True,
           cmap='viridis_r')


# In[24]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=["OTX1", "OTX2", "LMX1A", "EN1", "PITX2", "SIM2",
                  "HOXB-AS1", "HOTAIRM1", "HOXA2", "HOXB2", "GATA3", "GBX2",
                  "FGF8", "FGF17", "NKX2-8", "PAX8"],
           cmap='viridis_r')


# In[25]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.pl.umap(adata, color=['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', ],
           cmap='viridis_r')


# In[26]:


sc.set_figure_params(frameon=False, dpi=150, fontsize=8, figsize=(4, 3))
sc.tl.score_genes(adata, ["CALB1", "TH","SOX6"],
                      ctrl_size=50, gene_pool=None, n_bins=10, score_name='TH_score', random_state=0, copy=False, use_raw=True)
sc.pl.umap(adata, color='TH_score', cmap='viridis_r')
sc.pl.violin(
        adata,
        ["TH_score"],
        groupby="orig.ident",
        stripplot=False,  # remove the internal dots
        inner="box",  # adds a boxplot inside violins
    )


# In[30]:


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


# In[32]:


sc.pl.stacked_violin(
    adata, ['FOXA2', 'LMX1A', 'EN1', 'OTX2','SOX6','SHH','ASCL1','NEUROG2','DCX','NR4A2','TH','DDC','PITX3','CALB1','SLC18A2','SLC6A3'],
    groupby="seurat_clusters_pca_RNA_0.3", swap_axes=False, dendrogram=True
)


# In[ ]:





# In[33]:


sc.pl.umap(adata, color = ['orig.ident','seurat_clusters_pca_RNA_0.3'], legend_loc = 'on data',add_outline=True,  
                          legend_fontsize=5, legend_fontoutline=1,frameon=False,  size=25,)


# In[34]:


sc.pl.umap(adata, color = ['Phase'])


# In[ ]:




sc.pp.normalize_total(mid, target_sum = 1e4)
sc.pp.log1p(mid)sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
# In[35]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid_F.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[36]:


adata = predictions.to_adata()


# In[37]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_0.3', 'predicted_labels','majority_voting'])


# In[38]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA_0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA_0.3', use_as_prediction = 'majority_voting')


# In[ ]:





# In[39]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = '/Users/gclyu07/Desktop/AMS_analysis/celltypist_model_from_bimod_des_oc_F.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[40]:


adata = predictions.to_adata()
sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_0.3', 'predicted_labels','majority_voting'])


# In[41]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA_0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_pca_RNA_0.3', use_as_prediction = 'majority_voting')


# In[ ]:





# # 2. Label transferring proportion

# In[3]:


AMS <- readRDS('~/Desktop/AMS_analysis/0725_s2_AMS_ABV0.3_FILTERED.rds')


# In[7]:


cell_class <- unique(AMS$Cell_type_SL_abv0.3)
calculate_percentages <- function(x) {
    # Create a table with desired cell types, even if they are missing
  result_table <- table(factor(x, levels = cell_class))
  
  # Calculate proportions for specific cell types, including zeros
  cell_counts <- result_table[cell_class]
  missing_celltypes <- setdiff(cell_class, names(cell_counts))
  zero_counts <- rep(0, length(missing_celltypes))
  names(zero_counts) <- missing_celltypes
  final_counts <- c(cell_counts, zero_counts)
  
  prop.table(final_counts) * 100
}

percentages <- tapply(AMS$Cell_type_SL_abv0.3, AMS$seurat_clusters_pca_RNA_0.3, calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[8]:


labels


# In[9]:


table(AMS$orig.ident, AMS$seurat_clusters_pca_RNA_0.3)


# In[10]:


options(repr.plot.width=7, repr.plot.height=5)
DimPlot(AMS, group.by = 'seurat_clusters_pca_RNA_0.3', reduction = 'umap_pca_RNA_0.3', label = T)


# In[11]:


DimPlot(AMS, group.by = 'orig.ident', reduction = 'umap_pca_RNA_0.3', label = T)


# In[24]:


ml_oc_des <- c(
    "0_Rgl1", 
    "1_ProgFP.S",
    "2_NbM",
    "3_ProgFP.Development",
    "4_NbM",
    "5_ProgFP.Development",
    "6_ProgFP.Development",
    "7_ProgFP.G2M",
    "8_DA.Neuron projection",
    "9_DA.Synaptic assembly",
    "10_NProg",
    "11_ProgFP.G2M",
    "12_ProgFP.G2M",
    "13_Rgl2",
    "14_GabaNb",
    "15_NProg",
    "16_Rgl3", 
    "17_GabaNb",
    "18_DA.Synaptic assembly",
    "19_Rgl3",
    "20_ProgFP.S",
    "21_Rgl3",
    "22_NbM",
    "23_NProg"
)


# In[25]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)


# In[26]:


unique(ml_oc_des)


# In[27]:


length(unique(ml_oc_des))


# In[28]:


library(viridis)
cellType_merged_colors <-  viridis(11)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development", "NProg",
                                   "NbM","GabaNb", 
                                   "DA.Neuron projection","DA.Synaptic assembly")


# In[29]:


AMS$bimod_oc_des <- ml_oc_des[AMS@meta.data$seurat_clusters_pca_RNA_0.3]
AMS$bimod_oc_des <- factor(x = AMS$bimod_oc_des, levels = names(cellType_merged_colors))


# In[30]:


table(AMS$orig.ident,AMS$bimod_oc_des)


# In[31]:


options(repr.plot.width=10, repr.plot.height=7)
DimPlot(AMS, group.by = 'bimod_oc_des', ncol = 1, reduction = 'umap_pca_RNA', label = T)


# In[32]:


options(repr.plot.width=10, repr.plot.height=7)
DimPlot(AMS, group.by = 'bimod_oc_des', ncol = 1, reduction = 'umap_pca_RNA_0.3', label = T)


# In[39]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$orig.ident, AMS@meta.data$bimod_oc_des))
mat <- mat[,c("Rgl1",  "Rgl2", "Rgl3",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development", "NProg",
                                   "NbM","GabaNb", 
                                   "DA.Neuron projection","DA.Synaptic assembly")]
circos.clear()
par(cex = 1)
chordDiagram(mat, big.gap = 20, small.gap = 2, order = c(colnames(mat), rownames(mat)), grid.col = orig.ident_colors)  


# In[40]:


options(repr.plot.width=16, repr.plot.height=20)
VlnPlot(AMS, slot = "data", group.by='bimod_oc_des', features  =c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3"),
    pt.size = 0.2, ncol = 3)


# ### BI_MOD_OC

# In[42]:


ml_oc <- c(
    "0_Rgl1", 
    "1_ProgFP",
    "2_NbM",
    "3_ProgFP",
    "4_NbM",
    "5_ProgFP",
    "6_ProgFP",
    "7_ProgFP",
    "8_DA",
    "9_DA",
    "10_NProg",
    "11_ProgFP",
    "12_ProgFP",
    "13_Rgl2",
    "14_GabaNb",
    "15_NProg",
    "16_Rgl3", 
    "17_GabaNb",
    "18_DA",
    "19_Rgl3",
    "20_ProgFP",
    "21_Rgl3",
    "22_NbM",
    "23_NProg"
)


# In[43]:


ml_oc <- gsub(".*_", "", ml_oc)


# In[44]:


unique(ml_oc)


# In[45]:


length(unique(ml_oc))


# In[46]:


cellType_merged_colors <-  viridis(8)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3",
                                   "ProgFP", "NProg",
                                   "NbM","GabaNb", 
                                   "DA")


# In[47]:


AMS$bimod_oc <- ml_oc[AMS@meta.data$seurat_clusters_pca_RNA_0.3]
AMS$bimod_oc <- factor(x = AMS$bimod_oc, levels = names(cellType_merged_colors))


# In[48]:


saveRDS(AMS,"/Users/gclyu07/Desktop/AMS_analysis/0725_s3_AMS_annotated_v2.rds")


# # 3. MetaNeighbor

# In[3]:


mid <- readRDS('~/Desktop/XMAS_analysis/SL_midbrain_2.rds')


# In[4]:


AMS <- readRDS("/Users/gclyu07/Desktop/AMS_analysis/0725_s3_AMS_annotated_v2.rds")


# In[5]:


library(MetaNeighbor)
library(SummarizedExperiment)


# ## 3.1 Label transfer

# In[6]:


AMS@meta.data$dataset <- 'New_condition'
mid$dataset <- "Sten"


# In[7]:


Idents(mid)=mid@meta.data$LRprediction_labels
Idents(AMS)=AMS@meta.data$Cell_type_SL_abv0.3


# In[8]:


ob.list <- list(mid, AMS)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)
mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)
DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)

table(Idents(object = mda.integrated),mda.integrated@meta.data$dataset)
var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[9]:


Study_ID = rep(c('1', '2'), c(ncol(mid),ncol(AMS)))
Celltype = c(as.character(mid@meta.data$LRprediction_labels),as.character(AMS$Cell_type_SL_abv0.3))


# In[10]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[11]:


length(unique(mid@meta.data$LRprediction_labels))
length(unique(AMS$Cell_type_SL_abv0.3))


# In[12]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 18), rep('Experiment',15)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=12, repr.plot.height=12)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "average",
         annotation_row = ann_row,
         xtick  = FALSE)


# ### Neuronal reference

# In[13]:


mid_neu <- subset(mid, LRprediction_labels %in% c('DA','DA0','Gaba','GabaNb','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3'))


# In[14]:


Idents(mid_neu)=mid_neu@meta.data$LRprediction_labels
Idents(AMS)=AMS@meta.data$bimod_oc_des


# In[15]:


ob.list <- list(mid_neu, AMS)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)
mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)
DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)

table(Idents(object = mda.integrated),mda.integrated@meta.data$dataset)
var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[16]:


saveRDS(mda.integrated,"/Users/gclyu07/Desktop/AMS_analysis/outputs/0726_mid_neu_AMS_mda.integrated.rds")


# ______

# In[17]:


Study_ID = rep(c('1', '2'), c(ncol(mid_neu),ncol(AMS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(AMS$seurat_clusters_pca_RNA_0.3))


# In[18]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[19]:


length(unique(mid_neu@meta.data$LRprediction_labels))
length(unique(AMS$seurat_clusters_pca_RNA_0.3))


# In[20]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('New Condition',24)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=12, repr.plot.height=12)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "complete",
         annotation_row = ann_row,
         xtick  = FALSE)


# ______

# In[21]:


Study_ID = rep(c('1', '2'), c(ncol(mid_neu),ncol(AMS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(AMS$bimod_oc_des))


# In[22]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[23]:


length(unique(mid_neu@meta.data$LRprediction_labels))
length(unique(AMS$bimod_oc_des))


# In[24]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('New Condition',11)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=12, repr.plot.height=12)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "complete",
         annotation_row = ann_row,
         xtick  = FALSE)


# ______

# In[25]:


Study_ID = rep(c('1', '2'), c(ncol(mid_neu),ncol(AMS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(AMS$bimod_oc))


# In[26]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[27]:


length(unique(mid_neu@meta.data$LRprediction_labels))
length(unique(AMS$bimod_oc))


# In[28]:


table(AMS$bimod_oc)


# In[29]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('New Condition',8)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=12, repr.plot.height=12)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=3,cutree_cols=3,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "complete",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[ ]:





# # 4. Dotplots

# ## da marker genes

# In[121]:


options(repr.plot.width=10, repr.plot.height=6)
plot <- DotPlot(AMS, assay = "RNA", dot.min = 0.03, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "NEUROG2", "DCX"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[120]:


options(repr.plot.width=9, repr.plot.height=4)
plot <- DotPlot(AMS, assay = "RNA", group.by = 'bimod_oc', dot.min = 0.01, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6","KCNJ6","CALB1"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# ### cluster marker genes

# In[122]:


DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "seurat_clusters_pca_RNA_0.3"

AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[132]:


AMS$percent.rb <- PercentageFeatureSet(AMS, pattern = '^RP[SL]')
AMS$rbRatio <- AMS$percent.rb/100


# In[145]:


options(repr.plot.width=7.5, repr.plot.height=5)
VlnPlot(AMS, features = 'rbRatio', group.by = 'orig.ident')
options(repr.plot.width=20, repr.plot.height=5)
VlnPlot(AMS, features = 'rbRatio', group.by = 'seurat_clusters_pca_RNA_0.3')
VlnPlot(AMS, features = 'mitoRatio', group.by = 'seurat_clusters_pca_RNA_0.3')


# In[124]:


AMS.markers_best %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5


# In[162]:


filtered_genes <- top7[!top7$gene %like% '^RP[SL]',]


# In[166]:


options(repr.plot.width=50, repr.plot.height=10)
plot <- DotPlot(AMS, assay = "RNA", group.by = 'seurat_clusters_pca_RNA_0.3', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = unique(filtered_genes$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# ## cell type

# In[178]:


DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "bimod_oc"

AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[179]:


top <- AMS.markers_best %>% filter(!grepl("^RP", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[180]:


options(repr.plot.width=20, repr.plot.height=5)
plot <- DotPlot(AMS, assay = "RNA", group.by = 'bimod_oc', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "New condition", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# _____

# In[181]:


DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "bimod_oc_des"

AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[182]:


top <- AMS.markers_best %>% filter(!grepl("^RP", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[185]:


options(repr.plot.width=25, repr.plot.height=5)
plot <- DotPlot(AMS, assay = "RNA", group.by = 'bimod_oc_des', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 9, 
                features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "New condition", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# ______

# ## Volcano plot

# In[186]:


library(scRNAtoolVis)


# In[187]:


library(colorspace)
high_contrast_colors <- rainbow_hcl(11)


# In[195]:


deg <- FindAllMarkers(AMS, min.pct = 0.1, logfc.threshold = 0.5, only.pos = FALSE)


# In[229]:


options(repr.plot.width=19, repr.plot.height=10)
Idents(AMS) <- "bimod_oc_des"
jjVolcano(diffData = deg %>% filter(!grepl("^RP", gene)),
          tile.col = high_contrast_colors,
          pvalue.cutoff = 0.01) + NoLegend()


# In[250]:


options(repr.plot.width=16, repr.plot.height=8)
markerVocalno(markers = deg %>% filter(!grepl("^RP", gene)),labelCol = ggsci::pal_npg()(11)) + NoLegend()


# In[219]:


diff <- deg %>% filter(!grepl("^RP", gene)) %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[283]:


options(repr.plot.width=20, repr.plot.height=20)

jjDotPlot(object = AMS,
          gene = unique(diff$gene),
          id = 'bimod_oc_des',
          split.by = 'orig.ident',
          dot.col = c("#5CB85C", "#F0AD4E", "#D9534F"))


# In[230]:


options(repr.plot.width=20, repr.plot.height=10)
jjDotPlot(object = AMS,
          gene = unique(diff$gene),
          id = 'bimod_oc_des',
          dot.col = c("#5CB85C", "#D9534F"))


# In[234]:


options(repr.plot.width=20, repr.plot.height=10)
jjDotPlot(object = AMS,
          gene = unique(diff$gene),
          id = 'bimod_oc_des',
          ytree = F)


# In[279]:


options(repr.plot.width=24, repr.plot.height=5)
p1 <- DimPlot(AMS, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA", seed=seed, pt.size=0.5, label=TRUE)
p2 <- DimPlot(AMS, group.by="bimod_oc_des", reduction="umap_pca_RNA", seed=seed, pt.size=0.5, label=FALSE)
p3 <- DimPlot(AMS, group.by="orig.ident", reduction="umap_pca_RNA", seed=seed, pt.size=0.5, label=FALSE, cols = orig.ident_colors)
p1+p2+p3


# In[278]:


options(repr.plot.width=24, repr.plot.height=5)
p1 <- DimPlot(AMS, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=0.5, label=TRUE)
p2 <- DimPlot(AMS, group.by="bimod_oc_des", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=0.5, label=FALSE)
p3 <- DimPlot(AMS, group.by="orig.ident", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=0.5, label=FALSE, cols = orig.ident_colors)
p1+p2+p3


# # 5. ClusterProfiler

# In[42]:


AMS <- readRDS('~/Desktop/AMS_analysis/0725_s3_AMS_annotated_v2.rds')


# In[44]:


Idents(AMS) <- 'bimod_oc_des'
for (ct in unique(AMS$bimod_oc_des)) {
    
cluster.markers <- FindMarkers(AMS, ident.1 = ct,min.pct = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP[SL]", rownames(cluster.markers)))
    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  
options(repr.plot.width=7.2, repr.plot.height=3)
print(ego.bp %>% filter(p.adjust < 0.03) %>%
ggplot(showCategory = 6,
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
  scale_y_discrete(labels = label_wrap(45)) +
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}


# In[45]:


Idents(AMS) <- 'bimod_oc'
for (ct in unique(AMS$bimod_oc)) {
    
cluster.markers <- FindMarkers(AMS, ident.1 = ct,min.pct = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP[SL]", rownames(cluster.markers)))
    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  
options(repr.plot.width=7.2, repr.plot.height=3)
print(ego.bp %>% filter(p.adjust < 0.03) %>%
ggplot(showCategory = 6,
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
  scale_y_discrete(labels = label_wrap(45)) +
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}

