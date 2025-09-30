#!/usr/bin/env python
# coding: utf-8

# In[1]:


remove.packages('Seurat')


# In[2]:


remove.packages('SeuratObject')


# In[5]:


install.packages('~/Downloads/seurat-5.2.0.tar.gz')
install.packages('~/Downloads/seurat-object-5.1.0.tar.gz')


# In[6]:


install.packages('Seurat')


# In[7]:


install.packages('SeuratObject')


# # 0. Preprocessing

# In[1]:


XMAS <- readRDS('~/Desktop/XMAS_analysis/outputs/s06_MS_XMAS_relabeled.rds')


# In[8]:


XMAS <- readRDS('~/Desktop/XMAS_analysis/outputs/s03_XMAS_filtered_velocyto_RNA_BiModOC_labeled_with_unspliced.rds')


# In[9]:


XMAS


# In[10]:


library(scCustomize)


# In[11]:


library(viridis)
library(scRNAtoolVis)


# In[34]:


options(repr.plot.width=15, repr.plot.height=10)
FeaturePlot(XMAS, features = c('SULF1','SULF2','NDST3','NDST4'), reduction = 'umap', n=2)


# In[105]:


options(repr.plot.width=17.5, repr.plot.height=4)
plot <- DotPlot(XMAS, features = c("SULF1", "SULF2", "HS6ST1", "HS6ST2", "HS6ST3", "NDST1", "NDST2", "NDST3", "NDST4", "HS3ST1", "HS3ST2", "HS3ST4", "HS3ST5", "HS3ST6", "HS2ST1", "EXT1", "EXT2", "EXTL3", "XYLT1", "XYLT2", "GLCE"), 
        group.by   = 'bimod_oc',dot.min = 0.01, scale.by = 'size', scale = TRUE, dot.scale = 10,)
plot + theme(axis.text.x = element_text(angle = 90),  axis.text = element_text(size = 13)) + labs(x = "", y= "") +  scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[72]:


options(repr.plot.width=15, repr.plot.height=16)
VlnPlot(XMAS, features = c("SULF1", "SULF2", "HS6ST1", "HS6ST2", "HS6ST3", "NDST1", "NDST2", "NDST3", "NDST4", "HS3ST1", "HS3ST2", "HS3ST4", "HS3ST5", "HS3ST6", "HS2ST1", "EXT1", "EXT2", "EXTL3", "XYLT1", "XYLT2", "GLCE"), 
        group.by   = 'bimod_oc', n=4, alpha=0)


# In[ ]:


options(repr.plot.width=13, repr.plot.height=5)
plot <- DotPlot(sec.neu, assay = "RNA", group.by = 'celltype.split', dot.min = 0.01, scale.by= "size", scale=TRUE, dot.scale = 10,
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6","KCNJ6","CALB1","NEUROD1", "DCX","SYP","NRXN2","SREBF1"))
plot + theme(axis.text.x = element_text(angle = 90),  axis.text = element_text(size = 13)) +
labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[93]:


library(ggplot2)


# In[101]:


options(repr.plot.width=15, repr.plot.height=4.5)
plot <- DotPlot(XMAS, features = c("SDC1", "SDC2", "SDC3", "SDC4", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6", "HSPG2", "AGRN", "CD44", "COL18A1", "TBFBR3", "NRP1", "SRGN"), 
        dot.min = 0.01, scale.by = 'size', scale = TRUE, dot.scale = 10, group.by   = 'bimod_oc')
plot + theme(axis.text.x = element_text(angle = 90),  axis.text = element_text(size = 13)) + labs(x = "", y= "") +  scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[73]:


options(repr.plot.width=15, repr.plot.height=9)
VlnPlot(XMAS, features = c("SDC1", "SDC2", "SDC3", "SDC4", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6", "HSPG2", "AGRN", "CD44", "COL18A1", "TBFBR3", "NRP1", "SRGN"), 
        group.by   = 'bimod_oc',pt.size = NULL, alpha=0)


# In[76]:


options(repr.plot.width=15, repr.plot.height=6)
VlnPlot(XMAS, features = c("XYLT1", "XYLT2", "B4GALT7", "B3GALT6", "B3GAT3", "CSGalNAc-T1", "CSGalNAc-T2", "CHPF2"),
        group.by   = 'bimod_oc',pt.size = NULL,
  alpha = 0, n=4)


# In[78]:


options(repr.plot.width=15, repr.plot.height=3)
VlnPlot(XMAS, features = c("CHST11", "CHST3", "CHST15", "UST"),
        group.by   = 'bimod_oc',pt.size = NULL,
  alpha = 0, n=4)


# In[80]:


options(repr.plot.width=15, repr.plot.height=6)
VlnPlot(XMAS, features = c("DSE", "CHST14", "CHST12", "CHST3", "UST", "CHST15"),
        group.by   = 'bimod_oc',pt.size = NULL,
  alpha = 0, n=4)


# In[82]:


options(repr.plot.width=15, repr.plot.height=12)
VlnPlot(XMAS, features = c("CHST1", "CHST6", "B4GALT4", "B3GNT7", "B3GNT1", "B3GNT2", "B4GALT1", "B4GALT2", "B4GALT3", "ST3GAL1", "ST3GAL2", "ST3GAL3", "FUT8", "CHST2", "CHST4"),
        group.by   = 'bimod_oc',pt.size = NULL,
  alpha = 0, n=4)


# In[83]:


options(repr.plot.width=15, repr.plot.height=12)
VlnPlot(XMAS, features = c("IDS", "IDUA", "ARSB", "HYAL1", "HPSE", "HPSE2", "SGSH", "NAGLU", "HGSNAT", "GUSB", "GNS", "GALNS", "GLB1", "HEXA", "KERA"),
        group.by   = 'bimod_oc',pt.size = NULL,
  alpha = 0, n=4)


# In[ ]:





# In[22]:


VlnPlot(XMAS, features = c('SULF1','SULF2','NDST3','NDST4'), group.by   = 'orig.ident', n=2)


# In[14]:


XMAS


# In[13]:


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
    "15_ProgFP.Development",
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


# In[14]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)
unique(ml_oc_des)


# In[41]:


cellType_merged_colors <-  viridis(15)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3.ECM organization", "Rgl3.Purine metabolism", "Rgl3.Fate commitment",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development","NbM", "GabaNb", "Gaba",
                                   "DA.Neuron projection", "DA.Synaptic modulation","DA.Synaptic assembly","DA.Neurotransmitter release")
show_col(cellType_merged_colors)


# In[42]:


XMAS$bimod_oc_des <- ml_oc_des[XMAS@meta.data$seurat_clusters_BiMod_OC]
XMAS$bimod_oc_des <- factor(x = XMAS$bimod_oc_des, levels = names(cellType_merged_colors))


# In[43]:


options(repr.plot.width=15, repr.plot.height=10)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap", seed=seed, pt.size=2,label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_ATAC", seed=seed, pt.size=2, label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_BiMod", seed=seed, pt.size=2, label=TRUE)


# In[ ]:





# # 1. Celltypist loading

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


# # 2.  Transfer cell type labels from the STEN dataset

# In[15]:


sc.pp.normalize_total(mid, target_sum = 1e4)
sc.pp.log1p(mid)


# In[16]:


sc.pp.normalize_total(xmas, target_sum = 1e4)
sc.pp.log1p(xmas)


# ## 2.1 (Suggested) Feature selection for the first dataset

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


# ## 2.2 Model training and label transfer

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


# # 3. Mannual Annotation

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


# In[26]:


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
    "15_ProgFP",
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


# In[27]:


ml_oc <- gsub(".*_", "", ml_oc)


# In[28]:


unique(ml_oc)


# In[29]:


cellType_merged_colors <- c("#FF5733", "#33FF57", "#5733FF", "#FF33A6", "#33A6FF", "#A6FF33", "#A633FF", 
            "#FFC433", "#C433FF", "#FF334C", "#334CFF", "#FF3333")
names(cellType_merged_colors) <- c("DA",  "NbM", "Gaba", "GabaNb", "ProgFP", 
                                   "Rgl1", "Rgl2", "Rgl3")
show_col(cellType_merged_colors)


# In[30]:


XMAS$bimod_oc <- ml_oc[XMAS@meta.data$seurat_clusters_BiMod_OC]
XMAS$bimod_oc <- factor(x = XMAS$bimod_oc, levels = names(cellType_merged_colors))


# In[31]:


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


# In[37]:


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
    "15_ProgFP.Development", # changed from NbM
    "16_GabaNb", 
    "17_DA.Postsynaptic transmission",
    "18_DA.Postsynaptic transmission",
    "19_Rgl2",
    "20_Gaba",
    "21_ProgFP.S", 
    "22_DA.Neurotransmitter release",
    "23_NbM",
    "24_Rgl3.Purine metabolism",
    "25_Rgl1", 
    "26_ProgFP.G2M")


# In[38]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)
unique(ml_oc_des)


# In[39]:


cellType_merged_colors <-  viridis(15)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3.ECM organization", "Rgl3.Purine metabolism", "Rgl3.Fate commitment",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development","NbM", "GabaNb", "Gaba",
                                   "DA.Neuron projection", "DA.Synaptic development","DA.Postsynaptic transmission","DA.Neurotransmitter release")
show_col(cellType_merged_colors)


# In[40]:


XMAS$bimod_oc_des <- ml_oc_des[XMAS@meta.data$seurat_clusters_BiMod_OC]
XMAS$bimod_oc_des <- factor(x = XMAS$bimod_oc_des, levels = names(cellType_merged_colors))


# In[41]:


options(repr.plot.width=15, repr.plot.height=10)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap", seed=seed, pt.size=2,label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_ATAC", seed=seed, pt.size=2, label=TRUE)
DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap_BiMod", seed=seed, pt.size=2, label=TRUE)


# # 4. MetaNeighbor

# ## 4.1 Sten reference

# In[46]:


mid <- readRDS('~/Desktop/XMAS_analysis/SL_midbrain_2.rds')
mid$dataset <- "Sten"
Idents(mid)=mid@meta.data$LRprediction_labels
mid_neu <- subset(mid, LRprediction_labels %in% c('DA','DA0','Gaba','GabaNb','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3'))
XMAS@meta.data$dataset <- 'EA_proj'


# In[97]:


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


# In[43]:


mda.integrated <- readRDS('~/Desktop/XMAS_analysis/mda.integrated_XMAS_NEU_ABV031.rds')
table(Idents(object = mda.integrated),mda.integrated@meta.data$dataset)
var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[53]:


Study_ID = rep(c('1','2' ), c(ncol(mid_neu),ncol(XMAS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS$bimod_oc))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[56]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(Dataset=c(rep("Reference", 11), rep('Experiment',8)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=10, repr.plot.height=10)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=4,cutree_cols=4,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "complete",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[50]:


Study_ID = rep(c('1','2' ), c(ncol(mid_neu),ncol(XMAS)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS$bimod_oc_des))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[51]:


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


# ## 4.2 La Manno Ref.

# In[57]:


mid <- readRDS(file = '~/Desktop/XMAS_analysis/LaMannoEmbryo.rds')


# In[58]:


mid$dataset <- "La Manno"
mid <- as.Seurat(mid, data=NULL)

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


# In[60]:


library(Seurat)
library(SeuratDisk)


# In[61]:


SaveH5Seurat(XMAS, filename = ("~/Desktop/XMAS_analysis/outputs/s03_XMAS_filtered_velocyto_RNA_BiModOC_labeled_with_unspliced.h5Seurat"))
SeuratDisk::Convert("~/Desktop/XMAS_analysis/outputs/s03_XMAS_filtered_velocyto_RNA_BiModOC_labeled_with_unspliced.h5Seurat", dest = "h5ad")


# In[42]:


saveRDS(XMAS, '~/Desktop/XMAS_analysis/outputs/s03_XMAS_filtered_velocyto_RNA_BiModOC_labeled_with_unspliced.rds')


# # 5.Plots for publication

# ## 5.1 Proportion

# In[9]:


table_samples_by_cell_type <- XMAS@meta.data %>%
  dplyr::group_by(orig.ident, bimod_oc_des) %>%
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(bimod_oc_des, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))


# In[66]:


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


# In[69]:


options(repr.plot.width=8, repr.plot.height=7)

temp_labels <- seurat@meta.data %>%
  group_by(orig.ident) %>%
  tally()

table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(seurat@meta.data$orig.ident))) %>%
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

  scale_fill_manual(name = 'Cell type', values = sample(my36colors,15)) +
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


# In[110]:


table_samples_by_cell_type <- XMAS@meta.data %>%
  dplyr::group_by(orig.ident, bimod_oc) %>%
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(bimod_oc, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))


# In[111]:


my36colors <-c('#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


# In[114]:


options(repr.plot.width=8, repr.plot.height=7)

temp_labels <- XMAS@meta.data %>%
  group_by(orig.ident) %>%
  tally()

table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(XMAS@meta.data$orig.ident))) %>%
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

  scale_fill_manual(name = 'Cell type', values = sample(my36colors,15)) +
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


# In[ ]:





# In[149]:


table_samples_by_cell_type <- XMAS@meta.data %>%
  dplyr::group_by(orig.ident, bimod_oc) %>%
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(bimod_oc, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))


# In[178]:


my36colors <-c('#E63863', '#C1E6F3',  '#E0D4CA', '#5F3D69', '#C5DEBA',
                '#3A6963', '#DCC1DD', '#CCE0F5')


# In[182]:


options(repr.plot.width=6.5, repr.plot.height=7)

temp_labels <- XMAS@meta.data %>%
  group_by(orig.ident) %>%
  tally()

table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(XMAS@meta.data$orig.ident))) %>%
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

  scale_fill_manual(name = 'Cell type', values = my36colors) +
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


# In[27]:


colors <- c(
  rgb(0.651, 0.808, 0.890),  # Rgl1
  rgb(0.121, 0.466, 0.705),  # Rgl2
  rgb(0.698, 0.874, 0.541),  # Rgl3 ECM organization
  rgb(0.200, 0.627, 0.173),  # Rgl3 Purine metabolism
  rgb(0.984, 0.603, 0.600),  # Rgl3 Fate commitment
  rgb(0.890, 0.102, 0.110),  # ProgFP S
  rgb(0.992, 0.749, 0.435),  # ProgFP G2M
  rgb(1.000, 0.498, 0.000),  # ProgFP Development
  rgb(0.792, 0.698, 0.839),  # NbM
  rgb(0.415, 0.239, 0.603),  # GabaNb
  rgb(1.000, 1.000, 0.600),  # Gaba
  rgb(0.694, 0.349, 0.157),  # DA Neuron projection
  rgb(0.600, 0.600, 0.600),  # DA Synaptic assembly
  rgb(0.737, 0.502, 0.741),  # DA Neurotransmitter transmission
  rgb(0.698, 0.874, 0.541)   # DA Synaptic modulation (reuse from earlier)
)


# In[28]:


options(repr.plot.width=6.5, repr.plot.height=7)

temp_labels <- XMAS@meta.data %>%
  group_by(orig.ident) %>%
  tally()

table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(XMAS@meta.data$orig.ident))) %>%
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

  scale_fill_manual(name = 'Cell type', values = colors) +
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


# ## 5.2 Dot plots

# ## da marker genes

# In[44]:


Idents(XMAS) <- 'bimod_oc_des'


# In[45]:


XMAS <- XMAS[,!XMAS$bimod_oc_des %in% c('Gaba','GabaNb')]


# In[46]:


FeaturePlot(XMAS, features =c("NEUROD1"), reduction="umap",pt.size=1, label=TRUE)


# In[64]:


options(repr.plot.width=10, repr.plot.height=6)
plot <- DotPlot(XMAS, assay = "RNA", dot.min = 0.05, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6","KCNJ6","CALB1","NEUROD1", "DCX","SYP","NRXN2"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[48]:


options(repr.plot.width=10, repr.plot.height=6)
plot <- DotPlot(XMAS, assay = "RNA", dot.min = 0.03, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("SYP", "SYN1", "SYN2", "BSN", "SYNAP25","PSD95", "STXBP1"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[66]:


options(repr.plot.width=9, repr.plot.height=4)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'bimod_oc', dot.min = 0.01, scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6","KCNJ6","CALB1","NEUROD1", "DCX","SYP","NRXN2","SREBF1"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[ ]:





# ### cluster marker genes

# In[77]:


DefaultAssay(XMAS) <- "RNA"
Idents(object = XMAS) <- "seurat_clusters_BiMod_0.3"
 
XMAS.markers <- FindAllMarkers(XMAS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
XMAS.markers_best <- XMAS.markers[XMAS.markers$p_val_adj < 0.05,]


# In[78]:


XMAS$percent.rb <- PercentageFeatureSet(XMAS, pattern = '^RP[SL]')
XMAS$rbRatio <- XMAS$percent.rb/100


# In[79]:


options(repr.plot.width=7.5, repr.plot.height=5)
VlnPlot(XMAS, features = 'rbRatio', group.by = 'orig.ident')
options(repr.plot.width=20, repr.plot.height=5)
VlnPlot(XMAS, features = 'rbRatio', group.by = 'seurat_clusters_BiMod_0.3')
VlnPlot(XMAS, features = 'mitoRatio', group.by = 'seurat_clusters_BiMod_0.3')


# In[80]:


XMAS.markers_best %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5


# In[81]:


filtered_genes <- top5[!top5$gene %like% '^RP[SL]',]


# In[83]:


options(repr.plot.width=50, repr.plot.height=10)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'seurat_clusters_pca_RNA', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, dot.min = 0.1,
                features = unique(filtered_genes$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# ## cell type

# In[84]:


DefaultAssay(XMAS) <- "RNA"
Idents(object = XMAS) <- "bimod_oc_des"

XMAS.markers <- FindAllMarkers(XMAS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
XMAS.markers_best <- XMAS.markers[XMAS.markers$p_val_adj < 0.05,]


# In[85]:


top <- XMAS.markers_best %>% filter(!grepl("^RP", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[89]:


options(repr.plot.width=25, repr.plot.height=7.5)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'bimod_oc_des', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 8, dot.min = 0.1,
                features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[90]:


top <- XMAS.markers_best %>% filter(!grepl("^RP", gene)) %>% filter(!grepl("\\.", gene)) %>% filter(!grepl("-", gene)) %>% filter(!grepl("^LINC", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[92]:


options(repr.plot.width=25, repr.plot.height=7.5)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'bimod_oc_des', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 8, 
                 dot.min = 0.05, features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# _____

# In[93]:


DefaultAssay(XMAS) <- "RNA"
Idents(object = XMAS) <- "bimod_oc"

XMAS.markers <- FindAllMarkers(XMAS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
XMAS.markers_best <- XMAS.markers[XMAS.markers$p_val_adj < 0.05,]


# In[94]:


top <- XMAS.markers_best %>% filter(!grepl("^RP", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[99]:


options(repr.plot.width=15, repr.plot.height=5)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'bimod_oc_des', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 5, 
                dot.min = 0.05, features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = " ", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[100]:


top <- XMAS.markers_best %>% filter(!grepl("^RP", gene)) %>% filter(!grepl("\\.", gene)) %>% filter(!grepl("-", gene)) %>% filter(!grepl("^LINC", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[103]:


options(repr.plot.width=12.5, repr.plot.height=5)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'bimod_oc_des', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 5, 
                dot.min = 0.05, features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = " ", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# ______

# ## Volcano plot

# In[106]:


library(scRNAtoolVis)

library(colorspace)
high_contrast_colors <- rainbow_hcl(15)
# In[108]:


Idents(object = XMAS) <- "bimod_oc_des"
deg <- FindAllMarkers(XMAS, min.pct = 0.1, logfc.threshold = 0.5, only.pos = FALSE)


# In[117]:


options(repr.plot.width=40, repr.plot.height=15)
Idents(XMAS) <- "bimod_oc_des"
jjVolcano(diffData = deg %>% filter(!grepl("^RP", gene)),
          tile.col = sample(my36colors),
          pvalue.cutoff = 0.01) + NoLegend()


# In[119]:


options(repr.plot.width=40, repr.plot.height=20)
Idents(XMAS) <- "bimod_oc_des"
jjVolcano(diffData = deg %>% filter(!grepl("^RP", gene)) %>% 
          filter(!grepl("\\.", gene)) %>% filter(!grepl("-", gene)) %>% filter(!grepl("^LINC", gene)),
          tile.col = sample(my36colors),
          pvalue.cutoff = 0.01) + NoLegend()


# In[124]:


options(repr.plot.width=30, repr.plot.height=15)
markerVocalno(markers = deg %>% filter(!grepl("^RP", gene)),labelCol = sample(my36colors)) + NoLegend()


# In[126]:


options(repr.plot.width=25, repr.plot.height=10)
markerVocalno(markers = deg %>% filter(!grepl("^RP", gene)) %>% 
          filter(!grepl("\\.", gene)) %>% filter(!grepl("-", gene)) %>% filter(!grepl("^LINC", gene)),
              labelCol = sample(my36colors)) + NoLegend()


# In[127]:


diff <- deg %>% filter(!grepl("^RP", gene)) %>% 
          filter(!grepl("\\.", gene)) %>% filter(!grepl("-", gene)) %>% filter(!grepl("^LINC", gene))%>% filter(avg_log2FC > 0) %>% 
          group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[128]:


options(repr.plot.width=20, repr.plot.height=20)

jjDotPlot(object = XMAS,
          gene = unique(diff$gene),
          id = 'bimod_oc_des',
          split.by = 'orig.ident',
          dot.col = c("#5CB85C", "#F0AD4E", "#D9534F", "#66CCCC", "#AA98A9", "#9933CC"))


# In[129]:


options(repr.plot.width=20, repr.plot.height=10)
jjDotPlot(object = XMAS,
          gene = unique(diff$gene),
          id = 'bimod_oc_des',
          dot.col = c("#5CB85C", "#D9534F"))


# In[130]:


options(repr.plot.width=20, repr.plot.height=10)
jjDotPlot(object = XMAS,
          gene = unique(diff$gene),
          id = 'bimod_oc_des',
          ytree = F)


# In[131]:


options(repr.plot.width=20, repr.plot.height=10)
jjDotPlot(object = XMAS,
          gene = unique(diff$gene),
          id = 'bimod_oc',
          ytree = F)


# In[133]:


options(repr.plot.width=24, repr.plot.height=5)
p1 <- DimPlot(XMAS, group.by="seurat_clusters_BiMod_0.3", reduction="umap", seed=seed, pt.size=0.5, label=TRUE)
p2 <- DimPlot(XMAS, group.by="bimod_oc_des", reduction="umap", seed=seed, pt.size=0.5, label=FALSE)
p3 <- DimPlot(XMAS, group.by="orig.ident", reduction="umap", seed=seed, pt.size=0.5, label=FALSE, cols = orig.ident_colors)
p1+p2+p3


# # 5. ClusterProfiler

# In[135]:


Idents(XMAS) <- 'bimod_oc_des'
for (ct in unique(XMAS$bimod_oc_des)) {
    
cluster.markers <- FindMarkers(XMAS, ident.1 = ct,min.pct = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP[SL]", rownames(cluster.markers)) & 
                                               !grepl("\\.", rownames(cluster.markers)) & !grepl("-", rownames(cluster.markers)) & 
                                               !grepl("^LINC", rownames(cluster.markers)))

    
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
  scale_y_discrete(labels = label_wrap(45)) + NoLegend() +
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}


# In[136]:


Idents(XMAS) <- 'bimod_oc'
for (ct in unique(XMAS$bimod_oc)) {
    
cluster.markers <- FindMarkers(XMAS, ident.1 = ct,min.pct = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP[SL]", rownames(cluster.markers)) & 
                                               !grepl("\\.", rownames(cluster.markers)) & !grepl("-", rownames(cluster.markers)) & 
                                               !grepl("^LINC", rownames(cluster.markers)))

    
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
  scale_y_discrete(labels = label_wrap(45)) + NoLegend() +
theme(
  axis.text.y = element_text(size = 13),
  plot.title = element_text(size = 15,face="bold"),
  axis.text.x = element_text(size = 10),
  legend.text = element_text(size = 5),
  panel.grid.major = element_blank(),
))}


# In[137]:


Idents(XMAS) <- 'bimod_oc_des'
for (ct in unique(XMAS$bimod_oc_des)) {
cluster.markers <- FindMarkers(XMAS, ident.1 = ct,min.pct = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP[SL]", rownames(cluster.markers)) & 
                                               !grepl("\\.", rownames(cluster.markers)) & !grepl("-", rownames(cluster.markers)) & 
                                               !grepl("^LINC", rownames(cluster.markers)))

    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  

options(repr.plot.width=10, repr.plot.height=3.5)        

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


# In[139]:


Idents(XMAS) <- 'bimod_oc'
for (ct in unique(XMAS$bimod_oc)) {
cluster.markers <- FindMarkers(XMAS, ident.1 = ct,min.pct = 0.25)  

#RP gene filtering
cluster.markers <- cluster.markers %>% filter(!grepl("^RP[SL]", rownames(cluster.markers)) & 
                                               !grepl("\\.", rownames(cluster.markers)) & !grepl("-", rownames(cluster.markers)) & 
                                               !grepl("^LINC", rownames(cluster.markers)))

    
up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])
gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)  

options(repr.plot.width=10, repr.plot.height=3.5)        

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


# # 6. Python illustration

# In[1]:


import scanpy as sc
import scvelo as scv
import celltypist
import time
import numpy as np
import seaborn as sns


# In[2]:


adata = sc.read_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/s03_XMAS_filtered_velocyto_RNA_BiModOC_labeled_with_unspliced.h5ad")


# In[3]:


adata


# In[4]:


adata.obs['seurat_clusters_BiMod_0.3'] = adata.obs['seurat_clusters_BiMod_0.3'].astype(str)


# In[3]:


alist = {
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
12:'DA Synaptic assembly',
13:'DA Neurotransmitter transmission',
14:'DA Synaptic modulation'}
adata.obs['bimod_oc_des'] = (
adata.obs['bimod_oc_des']
.astype('category')
.map(alist)
)


# In[10]:


alist = {
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
.map(alist)
)


# In[14]:


sc.pl.umap(adata, color = ['bimod_oc_des'],frameon=False)


# In[16]:


sc.pl.umap(adata, color = ['orig.ident','seurat_clusters_BiMod_0.3'],add_outline=True,  
                          legend_fontsize=6, legend_fontoutline=1,frameon=False,  size=25,)


# In[8]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['orig.ident'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', palette=["#F0AD4E", "#D9534F", "#428BCA", "#9933CC", "#66CCCC"])


# In[11]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['bimod_oc'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', 
           palette=['#E63863', '#C1E6F3',  '#E0D4CA', '#5F3D69', '#C5DEBA', '#3A6963', '#DCC1DD', '#CCE0F5'])


# In[34]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['bimod_oc_des'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', 
           palette=['#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59',])


# In[6]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['bimod_oc_des'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='')


# In[19]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['bimod_oc_des'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', palette='Paired')


# In[14]:


palette = sns.color_palette("Paired", 36) 


# In[18]:


palette = [
    (1.000, 1.000, 0.600),  # Rgl1
    (0.694, 0.349, 0.157),  # Rgl2
    (0.600, 0.600, 0.600),  # Rgl3 ECM organization
    (0.737, 0.502, 0.741),  # Rgl3 Purine metabolism
    (0.698, 0.874, 0.541),  # Rgl3 Fate commitment



    (0.415, 0.239, 0.603),  # ProgFP S
    (0.792, 0.698, 0.839),  # ProgFP G2M
    (1.000, 0.498, 0.000),  # ProgFP Development

    (0.992, 0.749, 0.435),  # NbM
    (0.984, 0.603, 0.600),  # Gaba
    (0.890, 0.102, 0.110),  # GabaNb
    (0.651, 0.808, 0.890),  # DA Neuron projection
    (0.121, 0.466, 0.705),  # DA Synaptic assembly
    (0.698, 0.874, 0.541),  # DA Neurotransmitter transmission
    (0.200, 0.627, 0.173),  # DA Synaptic modulation
]


# In[19]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['bimod_oc_des'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', palette=palette)


# In[17]:


palette


# In[52]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color = ['bimod_oc'],add_outline=True,  
                          legend_fontsize=10, legend_fontoutline=1,frameon=False,  size=50, title='', palette='Paired')


# In[ ]:





# In[ ]:





# # 7. Rerun the dataset with celltypist (?)

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

