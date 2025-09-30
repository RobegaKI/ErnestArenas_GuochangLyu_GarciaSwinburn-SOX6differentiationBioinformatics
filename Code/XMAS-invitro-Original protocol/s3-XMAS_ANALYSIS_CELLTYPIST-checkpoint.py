#!/usr/bin/env python
# coding: utf-8

# XMAS_ANALYSIS_05_CELLTYPIST 051424

# In[665]:


ml_oc_des <- c(
    "0_Rgl1", 
    "1_Rgl2",
    "2_ProgFP.Development",
    "3_Rgl3.ECM organization",
    "4_ProgFP.S",
    "5_ProgFP.G2M",
    "6_DA.Neurotransmitter release",
    "7_DA.Neurotransmitter release",
    "8_DA.Neuron projection",
    "9_ProgFP.G2M",
    "10_DA.Synaptic assembly",
    "11_DA.Neuron projection",
    "12_DA.Synaptic modulation",
    "13_DA.Synaptic modulation",
    "14_Rgl3.Fate commitment",
    "15_NbM",
    "16_GabaNb", 
    "17_DA.Synaptic assembly",
    "18_DA.Synaptic assembly",
    "19_Rgl2",
    "20_Gaba",
    "21_ProgFP.S", 
    "22_DA.Neurotransmitter release",
    "23_NbM",
    "24_Rgl3.Purine metabolism",
    "25_Rgl1", 
    "26_ProgFP.G2M")


# In[15]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)
unique(ml_oc_des)


# In[16]:


cellType_merged_colors <-  viridis(15)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3.ECM organization", "Rgl3.Purine metabolism", "Rgl3.Fate commitment",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development","NbM", "GabaNb", "Gaba",
                                   "DA.Neuron projection", "DA.Synaptic modulation","DA.Synaptic assembly","DA.Neurotransmitter release")
show_col(cellType_merged_colors)


# In[17]:


XMAS$bimod_oc_des <- ml_oc_des[XMAS@meta.data$seurat_clusters_BiMod_OC]
XMAS$bimod_oc_des <- factor(x = XMAS$bimod_oc_des, levels = names(cellType_merged_colors))


# In[94]:


options(repr.plot.width=15, repr.plot.height=10)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap", seed=seed, pt.size=2,label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_ATAC", seed=seed, pt.size=2, label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_BiMod", seed=seed, pt.size=2, label=TRUE)


# In[ ]:





# In[2]:


library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(sctransform)
library(dplyr)
library(stringr)
library(clustree)
library(data.table)

library(cowplot)
library(randomForest)
library(ggplot2)
library(caret)

library(biomaRt)

library(pheatmap)
library(grid)
library(gridExtra)
library(biovizBase)
library(ggrepel)
library(viridis)
library(reshape2)
library(hues)
library(clusterProfiler)
library(org.Hs.eg.db)

library(tidyverse)
library(ggnewscale)
library(circlize)
library(scales)
library(hrbrthemes)
library(patchwork)
library(cluster)
library(corrplot)
library(qgraph)
library(plyr)

library(SummarizedExperiment)
library(gplots)
library(RColorBrewer)


# In[1]:


XMAS <- readRDS('~/Desktop/XMAS_analysis/03_abv031_BiMerged_XMAS.rds')


# In[3]:


XMAS <- FindMultiModalNeighbors(
  object = XMAS,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:25),
  weighted.nn.name = "weighted.nn",
  modality.weight.name = c("RNA.weight.name", "peaks.weight.name"),
  snn.graph.name = "wsnn",
  prune.SNN = 0,
  verbose = TRUE
)


# In[4]:


XMAS <- FindClusters(XMAS, resolution = 3.3, verbose = FALSE, graph.name="wsnn")
XMAS$seurat_clusters_BiMod_OC <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# In[5]:


table(XMAS$seurat_clusters_BiMod_OC)


# In[10]:


table(XMAS$seurat_clusters_BiMod_OC)


# In[6]:


XMAS


# In[7]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(XMAS, group.by="seurat_clusters_BiMod_OC", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[12]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(XMAS, group.by="seurat_clusters_BiMod_OC", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[120]:


library(Seurat)
library(SeuratDisk)


# In[20]:


SaveH5Seurat(XMAS, filename = ("~/Desktop/XMAS_analysis/outputs/11_XMAS_filtered_velocyto_RNA_BiModOC.h5Seurat"))
SeuratDisk::Convert("~/Desktop/XMAS_analysis/outputs/11_XMAS_filtered_velocyto_RNA_BiModOC.h5Seurat", dest = "h5ad")


# In[ ]:





# In[ ]:





# In[1]:


import scanpy as sc
import scvelo as scv
import celltypist
import time
import numpy as np


# In[2]:


mid = sc.read_h5ad('Desktop/XMAS_analysis/outputs_no_d0/SL_midbrain_2.h5ad',)


# In[265]:


xmas = sc.read_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/11_XMAS_filtered_velocyto_RNA_BiModOC.h5ad")


# In[27]:


xmas_unspliced = sc.read_h5ad("Desktop/XMAS_analysis/outputs_no_d0/11_XMAS_filtered_velocyto_unspliced.h5ad")


# In[12]:


xmas


# In[266]:


list = {
0:'DA',
1:'DA0',
2:'NbML1',
3:'NbM',
4:'NProg',
5:'ProgFP',
6:'Rgl1',
7:'Rgl2',
8:'Rgl3',
}
xmas.obs['BiMod_merged_0.3'] = (
xmas.obs['BiMod_merged_0.3']
.map(list)
.astype('category')
)


# In[22]:


xmas_unspliced


# In[14]:


xmas.raw.X


# In[134]:


xmas.X


# In[59]:


xmas.layers["unspliced"]


# In[58]:


xmas


# In[36]:


mid


# ## ## Transfer cell type labels from the first dataset to the second dataset
# 

# In[15]:


sc.pp.normalize_total(mid, target_sum = 1e4)
sc.pp.log1p(mid)


# In[16]:


sc.pp.normalize_total(xmas, target_sum = 1e4)
sc.pp.log1p(xmas)


# ### ### (Suggested) Feature selection for the first dataset

# In[187]:


t_start = time.time()
model_fs = celltypist.train(mid[:,:], 'LRprediction_labels', n_jobs = 10, max_iter = 5, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[188]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -100, axis = 1)[:, -100:]


# In[189]:


gene_index = np.unique(gene_index)


# In[190]:


print(f"Number of genes selected: {len(gene_index)}")


# ### ### Model training and label transfer

# In[203]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid[:,gene_index], 'LRprediction_labels', check_expression = False, n_jobs = 10, max_iter = 150)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[204]:


model.write('/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid_3.pkl')


# In[205]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(xmas, model = '/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid_3.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[206]:


adata = predictions.to_adata()


# In[207]:


adata.obs.iloc[:, -4:]


# In[97]:


xmas.raw.X


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### ### (Suggested) Feature selection for the first dataset

# In[153]:


t_start = time.time()
model_fs = celltypist.train(mid[:,:], 'LRprediction_labels', n_jobs = 10, max_iter = 5, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[166]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -200, axis = 1)[:, -200:]


# In[167]:


gene_index = np.unique(gene_index)


# In[168]:


print(f"Number of genes selected: {len(gene_index)}")


# ### ### Model training and label transfer

# In[169]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid[:,gene_index], 'LRprediction_labels', check_expression = False, n_jobs = 10, max_iter = 100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[170]:


model.write('/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid.pkl')


# In[171]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(xmas, model = '/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[172]:


predictions.predicted_labels


# In[173]:


celltypist.dotplot(predictions, use_as_reference = 'Cell_type_SL_abv0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'Cell_type_SL_abv0.3', use_as_prediction = 'majority_voting')


# In[174]:


celltypist.dotplot(predictions, use_as_reference = 'BiMod_merged_0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'BiMod_merged_0.3', use_as_prediction = 'majority_voting')


# In[175]:


adata = predictions.to_adata()


# In[176]:


adata.obs.iloc[:, -4:]


# In[177]:


sc.pl.umap(xmas, color = ['BiMod_merged_0.3', 'predicted_labels','majority_voting'])


# In[ ]:





# ### ### (Suggested) Feature selection for the first dataset

# In[213]:


t_start = time.time()
model_fs = celltypist.train(mid[:,:], 'LRprediction_labels', n_jobs = 10, max_iter = 5, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[310]:


gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -130, axis = 1)[:, -130:]


# In[311]:


gene_index = np.unique(gene_index)


# In[312]:


print(f"Number of genes selected: {len(gene_index)}")


# ### ### Model training and label transfer

# In[313]:


# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(mid[:,gene_index], 'LRprediction_labels', check_expression = False, n_jobs = 10, max_iter = 150)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")


# In[314]:


model.write('/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid_F.pkl')


# In[315]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(xmas, model = '/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid_F.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[316]:


adata = predictions.to_adata()
sc.pl.umap(xmas, color = ['BiMod_merged_0.3', 'predicted_labels','majority_voting'])


# In[317]:


predictions.predicted_labels


# In[318]:


celltypist.dotplot(predictions, use_as_reference = 'Cell_type_SL_abv0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'Cell_type_SL_abv0.3', use_as_prediction = 'majority_voting')


# In[319]:


celltypist.dotplot(predictions, use_as_reference = 'BiMod_merged_0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'BiMod_merged_0.3', use_as_prediction = 'majority_voting')


# In[320]:


adata.obs.iloc[:, -4:]


# In[336]:


pd.crosstab(adata.obs['predicted_labels'],adata.obs['seurat_clusters_BiMod_OC'])


# In[337]:


pd.crosstab(adata.obs['majority_voting'],adata.obs['seurat_clusters_BiMod_OC'])


# In[334]:


pl.to_csv('/Users/gclyu07/Desktop/XMAS_analysis/outputs/celltypist_F_pl.csv')


# In[335]:


mv.to_csv('/Users/gclyu07/Desktop/XMAS_analysis/outputs/celltypist_F_mv.csv')


# In[ ]:





# In[20]:


cell_class <- unique(XMAS$Cell_type_SL)
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

percentages <- tapply(XMAS$Cell_type_SL, XMAS@meta.data[["seurat_clusters_BiMod_OC"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[21]:


labels


# In[25]:


write_csv(labels, "~/Desktop/XMAS_analysis/outputs/BiMod_OC.csv")


# In[15]:


table(XMAS$Cell_type_SL,XMAS$seurat_clusters_BiMod_OC)


# In[62]:


table(XMAS$seurat_clusters_BiMod_OC)


# In[12]:


XMAS


# In[98]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(XMAS, group.by="seurat_clusters_BiMod_OC", reduction="umap", seed=seed, pt.size=1, label=TRUE)
DimPlot(XMAS, group.by="orig.ident", reduction="umap", seed=seed, pt.size=1, label=TRUE)
DimPlot(XMAS, group.by="Phase", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[8]:


ml_oc<- c(
    "0_Rgl1", 
    "1_Rgl2",
    "2_ProgFP",
    "3_Rgl3",
    "4_ProgFP",
    "5_ProgFP",
    "6_DA",
    "7_DA",
    "8_DA",
    "9_ProgFP",
    "10_DA",
    "11_DA",
    "12_DA",
    "13_DA",
    "14_Rgl3",
    "15_NbM",
    "16_GabaNb", 
    "17_DA",
    "18_DA",
    "19_Rgl2",
    "20_Gaba",
    "21_ProgFP", 
    "22_DA",
    "23_NbM",
    "24_Rgl3",
    "25_Rgl1", 
    "26_ProgFP")


# In[9]:


ml_oc <- gsub(".*_", "", ml_oc)


# In[10]:


unique(ml_oc)


# In[11]:


cellType_merged_colors <- c("#FF5733", "#33FF57", "#5733FF", "#FF33A6", "#33A6FF", "#A6FF33", "#A633FF", 
            "#FFC433", "#C433FF", "#FF334C", "#334CFF", "#FF3333")
names(cellType_merged_colors) <- c("DA",  "NbM", "Gaba", "GabaNb", "ProgFP", 
                                   "Rgl1", "Rgl2", "Rgl3")
show_col(cellType_merged_colors)


# In[12]:


XMAS$bimod_oc <- ml_oc[XMAS@meta.data$seurat_clusters_BiMod_OC]
XMAS$bimod_oc <- factor(x = XMAS$bimod_oc, levels = names(cellType_merged_colors))


# In[87]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(XMAS, group.by="orig.ident", reduction="umap", seed=seed, pt.size=1, label=TRUE)

DimPlot(XMAS, group.by="bimod_oc", reduction="umap", seed=seed, pt.size=1, label=TRUE)
DimPlot(XMAS, group.by="bimod_oc", reduction="umap_ATAC", seed=seed, pt.size=1, label=TRUE)
DimPlot(XMAS, group.by="bimod_oc", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[45]:


DefaultAssay(XMAS) <- "RNA"
Idents(XMAS) <- "seurat_clusters_BiMod_OC"


# In[48]:


DefaultAssay(XMAS) <- "RNA"

XMAS.markers <- FindAllMarkers(XMAS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
XMAS.markers_best <- XMAS.markers[XMAS.markers$p_val_adj < 0.05,]
topDiffCluster <- XMAS.markers_best %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


# In[49]:


sample_df <- as.data.frame(lapply(split(topDiffCluster, topDiffCluster$cluster), function(x) c(x$gene,rep("None", 100-length(x$gene)))))
colnames(sample_df) <- levels(topDiffCluster$cluster)
sample_df


# In[51]:


write.csv(sample_df,"~/Desktop/XMAS_analysis/outputs/05_XMAS_Top100_HVGs_0.3_BimodOC.csv",row.names = FALSE)


# In[14]:


ml_oc_des <- c(
    "0_Rgl1", 
    "1_Rgl2",
    "2_ProgFP.Development",
    "3_Rgl3.ECM organization",
    "4_ProgFP.S",
    "5_ProgFP.G2M",
    "6_DA.Neurotransmitter release",
    "7_DA.Neurotransmitter release",
    "8_DA.Neuron projection",
    "9_ProgFP.G2M",
    "10_DA.Postsynaptic transmission",
    "11_DA.Neuron projection",
    "12_DA.Synaptic development",
    "13_DA.Synaptic development",
    "14_Rgl3.Fate commitment",
    "15_NbM",
    "16_GabaNb", 
    "17_DA.Postsynaptic transmission",
    "18_DA.Postsynaptic transmission",
    "19_Rgl2",
    "20_Gaba",
    "21_ProgFP.S", 
    "22_DA.Neurotransmitter release",
    "23_cNbM",
    "24_Rgl3.Purine metabolism",
    "25_Rgl1", 
    "26_ProgFP.G2M")


# In[15]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)
unique(ml_oc_des)

cellType_merged_colors <-  viridis(15)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3.ECM organization", "Rgl3.Purine metabolism", "Rgl3.Fate commitment",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development","NbM", "GabaNb", "Gaba",
                                   "DA.Neuron projection", "DA.Synaptic development","DA.Postsynaptic transmission","DA.Neurotransmitter release")
show_col(cellType_merged_colors)
# In[17]:


XMAS$bimod_oc_des <- ml_oc_des[XMAS@meta.data$seurat_clusters_BiMod_OC]
XMAS$bimod_oc_des <- factor(x = XMAS$bimod_oc_des, levels = names(cellType_merged_colors))


# In[94]:


options(repr.plot.width=15, repr.plot.height=10)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap", seed=seed, pt.size=2,label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_ATAC", seed=seed, pt.size=2, label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_BiMod", seed=seed, pt.size=2, label=TRUE)


# In[ ]:





# In[97]:


mid <- readRDS('~/Desktop/XMAS_analysis/SL_midbrain_2.rds')
mid$dataset <- "Sten"
Idents(mid)=mid@meta.data$LRprediction_labels
mid_neu <- subset(mid, LRprediction_labels %in% c('DA','DA0','Gaba','GabaNb','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3'))
XMAS@meta.data$dataset <- 'EA_proj'
ob.list <- list(mid_neu, XMAS)
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

Study_ID = rep(c('1','2' ), c(ncol(mid_neu),ncol(XMAS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS$bimod_oc))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',8)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=10, repr.plot.height=10)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "complete",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )
# In[118]:


Study_ID = rep(c('1','2' ), c(ncol(mid_neu),ncol(XMAS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS$bimod_oc_des))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[119]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(Dataset=c(rep("Reference", 11), rep('Experiment',15)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=10, repr.plot.height=10)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "complete",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[122]:


mid <- readRDS(file = '~/Desktop/XMAS_analysis/LaMannoEmbryo.rds')


# In[123]:


mid$dataset <- "La Manno"
mid <- as.Seurat(mid, data=NULL)


# In[124]:


ob.list <- list(mid, XMAS)

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


# In[143]:


Study_ID = rep(c('1','2' ), c(ncol(mid),ncol(XMAS)))
Celltype = c(as.character(mid@meta.data$Cell_type),as.character(XMAS$bimod_oc))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[144]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 26), rep('Experiment',8)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=10, repr.plot.height=10)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[145]:


Study_ID = rep(c('1','2' ), c(ncol(mid),ncol(XMAS)))
Celltype = c(as.character(mid@meta.data$Cell_type),as.character(XMAS$bimod_oc_des))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[157]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(Dataset=c(rep("Reference", 26), rep('Experiment',15)))
rownames(ann_row)=rownames(celltype_NV)

options(repr.plot.width=15, repr.plot.height=10)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         ) 


# In[19]:


saveRDS(XMAS,"~/Desktop/XMAS_analysis/outputs/05_XMAS_filtered_velocyto_RNA_BiModOC_labeled.rds")


# In[121]:


SaveH5Seurat(XMAS, filename = ("~/Desktop/XMAS_analysis/outputs/05_XMAS_filtered_velocyto_RNA_BiModOC_labeled.h5Seurat"))
SeuratDisk::Convert("~/Desktop/XMAS_analysis/outputs/05_XMAS_filtered_velocyto_RNA_BiModOC_labeled.h5Seurat", dest = "h5ad")


# In[ ]:





# In[ ]:





# In[ ]:





# In[477]:


import scanpy as sc
import scvelo as scv
import celltypist
import time
import numpy as np


# In[501]:


xmas = sc.read_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/11_XMAS_filtered_velocyto_RNA_BiModOC.h5ad")


# In[502]:


xmas.X


# ## ## Transfer cell type labels from the first dataset to the second dataset
# 

# In[492]:


sc.pp.normalize_total(xmas, target_sum = 1e4)
sc.pp.log1p(xmas)


# ### ### Model training and label transfer

# In[503]:


# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(xmas, model = '/Users/gclyu07/Desktop/XMAS_analysis/celltypist_model_from_mid_F.pkl',majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")


# In[504]:


adata = predictions.to_adata()
sc.pl.umap(xmas, color = ['BiMod_merged_0.3', 'predicted_labels','majority_voting'])


# In[497]:


adata = predictions.to_adata()
sc.pl.umap(xmas, color = ['BiMod_merged_0.3', 'predicted_labels','majority_voting'])


# In[498]:


predictions.predicted_labels


# In[505]:


adata


# In[506]:


celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_BiMod_OC', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters_BiMod_OC', use_as_prediction = 'majority_voting')


# In[319]:


celltypist.dotplot(predictions, use_as_reference = 'BiMod_merged_0.3', use_as_prediction = 'predicted_labels')
celltypist.dotplot(predictions, use_as_reference = 'BiMod_merged_0.3', use_as_prediction = 'majority_voting')


# In[320]:


adata.obs.iloc[:, -4:]


# In[336]:


pd.crosstab(adata.obs['predicted_labels'],adata.obs['seurat_clusters_BiMod_OC'])


# In[337]:


pd.crosstab(adata.obs['majority_voting'],adata.obs['seurat_clusters_BiMod_OC'])

