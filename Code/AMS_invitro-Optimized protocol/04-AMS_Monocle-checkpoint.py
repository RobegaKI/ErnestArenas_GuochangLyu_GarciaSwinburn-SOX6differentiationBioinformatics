#!/usr/bin/env python
# coding: utf-8

# # Monocle

# In[74]:


BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))


# In[9]:


devtools::install_github('cole-trapnell-lab/monocle3',force = TRUE)


# In[10]:


remotes::install_github('satijalab/seurat-wrappers')


# In[24]:


install.packages('Seurat')


# In[1]:


library(monocle3)
library(SeuratWrappers)


# In[ ]:





# https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed

# # Monocle analyzing (DA lineage)

# In[2]:


library(Seurat)
library(monocle3)


# In[3]:


sec <- readRDS('~/Desktop/AMS_analysis/0725_s3_AMS_annotated.rds')


# In[4]:


table(sec$orig.ident, sec$bimod_oc)


# ## Pseudotime arrangement (confirmed)

# In[549]:


ml_oc_des <- c(
    "0_Rgl1", 
    "1_ProgFP.G2M",
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
    "20_ProgFP.G2M",
    "21_Rgl3",
    "22_NbM",
    "23_NProg"
)


# In[550]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)


# In[551]:


unique(ml_oc_des)


# In[552]:


cellType_merged_colors <-  viridis(11)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development", "NProg",
                                   "NbM","GabaNb", 
                                   "DA.Neuron projection","DA.Synaptic assembly")


# In[553]:


sec$bimod_oc_des <- ml_oc_des[sec@meta.data$seurat_clusters_pca_RNA_0.3]
sec$bimod_oc_des <- factor(x = sec$bimod_oc_des, levels = names(cellType_merged_colors))


# In[554]:


cds <- subset(sec, subset = bimod_oc %in% c("DA","NbM","NProg","ProgFP") 
              & !(sec$bimod_oc_des %in% c('ProgFP.G2M', 'ProgFP.S','DA.Neuron projection'))
              & !(sec$bimod_oc == "ProgFP" & sec$orig.ident == 'New-D28'))


# In[555]:


cds


# In[556]:


options(repr.plot.width=18, repr.plot.height=4)
p1<-DimPlot(cds, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1,label=TRUE)
p2<-DimPlot(cds, group.by="bimod_oc_des", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1,label=FALSE)
p3<-DimPlot(cds, group.by="bimod_oc", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=TRUE)
p1+p2+p3


# In[557]:


cds <- as.cell_data_set(cds)


# In[558]:


cds <- preprocess_cds(cds)


# In[559]:


cds <- reduce_dimension(cds)


# In[560]:


options(repr.plot.width=5, repr.plot.height=3)
monocle3::plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
monocle3::plot_cells(cds, color_cells_by="bimod_oc", label_cell_groups=FALSE)


# In[561]:


cds <- cluster_cells(cds, resolution = 1e-3)


# In[562]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds,  color_cells_by = "bimod_oc_des", label_groups_by_cluster=FALSE)


# In[563]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster",
           label_cell_groups = TRUE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")


# In[564]:


rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name


# In[565]:


genes <- c("FOXA2","LMX1A","TH","NR4A2","EN1","SOX6")


# In[566]:


# patterns of expression
options(repr.plot.width=10, repr.plot.height=5)
plot_cells(cds,
           genes=genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


# In[567]:


cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)


# In[568]:


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) %in% c(5,1)]))


# In[569]:


cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "orig.ident") 


# In[570]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


# In[571]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "bimod_oc_des",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")


# In[572]:


pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)


# In[573]:


library(ggbeeswarm)


# In[574]:


library(RColorBrewer)
par(mar=c(3,4,2,2))
options(repr.plot.width=8, repr.plot.height=10)
display.brewer.all()


# In[575]:


options(repr.plot.width=8, repr.plot.height=3)
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = orig.ident, colour = bimod_oc)) + 
    geom_quasirandom(groupOnX = FALSE) + theme_classic() + 
    scale_color_brewer(palette = "Paired")  +
    xlab("Pseudotime") + ylab("Cell types") +  labs(colour = " ") 


# In[576]:


options(repr.plot.width=7, repr.plot.height=3)
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = bimod_oc, colour = orig.ident)) + 
    geom_quasirandom(groupOnX = FALSE) + theme_classic() + 
    scale_color_manual(values = c("#5CB85C", "#F0AD4E", "#D9534F")) +
    xlab("Pseudotime") + ylab("Cell types") +  labs(colour = " ") 


# In[577]:


head(pseudotime(cds), 10)


# In[578]:


cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
options(repr.plot.width=8, repr.plot.height=4)
ggplot(data.pseudo, aes(monocle3_pseudotime, orig.ident, fill = orig.ident)) + geom_boxplot()


# In[579]:


ggplot(data.pseudo, aes(monocle3_pseudotime, orig.ident, fill = bimod_oc)) + geom_boxplot()


# In[580]:


ggplot(data.pseudo, aes(monocle3_pseudotime, bimod_oc, fill = orig.ident)) + geom_boxplot()


# In[581]:


options(repr.plot.width=8, repr.plot.height=6)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters_pca_RNA, monocle3_pseudotime), fill = seurat_clusters_pca_RNA)) + geom_boxplot()


# In[582]:


cds <- estimate_size_factors(cds)


# In[583]:


options(repr.plot.width=4, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime") 


# In[584]:


options(repr.plot.width=3.5, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "bimod_oc") 


# In[585]:


options(repr.plot.width=3.5, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "orig.ident") 


# ____

# ## Subset alternative (w. cluster projection)

# In[212]:


cds <- subset(sec, subset = bimod_oc %in% c("DA","NbM","NProg","ProgFP") 
              & !(sec$bimod_oc_des %in% c('ProgFP.G2M', 'ProgFP.S'))
              & !(sec$bimod_oc == "ProgFP" & sec$orig.ident == 'New-D28'))


# In[213]:


cds


# In[214]:


options(repr.plot.width=18, repr.plot.height=4)
p1<-DimPlot(cds, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1,label=TRUE)
p2<-DimPlot(cds, group.by="bimod_oc_des", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1,label=FALSE)
p3<-DimPlot(cds, group.by="bimod_oc", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=TRUE)
p1+p2+p3


# In[215]:


cds <- as.cell_data_set(cds)


# In[216]:


cds <- preprocess_cds(cds)


# In[217]:


cds <- reduce_dimension(cds)


# In[218]:


options(repr.plot.width=5, repr.plot.height=3)
monocle3::plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
monocle3::plot_cells(cds, color_cells_by="bimod_oc", label_cell_groups=FALSE)


# In[219]:


cds <- cluster_cells(cds, resolution = 1e-3)


# In[220]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds,  color_cells_by = "bimod_oc_des", label_groups_by_cluster=FALSE)


# In[221]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster",
           label_cell_groups = TRUE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")


# In[222]:


rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name


# In[225]:


cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)


# In[453]:


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) %in% c(1,6,4,19)]))


# In[227]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


# In[228]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "bimod_oc_des",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")


# In[229]:


pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)


# In[230]:


options(repr.plot.width=8, repr.plot.height=3)
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = orig.ident, colour = bimod_oc)) + 
    geom_quasirandom(groupOnX = FALSE) + theme_classic() + 
    scale_color_brewer(palette = "Paired")  +
    xlab("Pseudotime") + ylab("Cell types") +  labs(colour = " ") 


# In[231]:


options(repr.plot.width=7, repr.plot.height=3)
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = bimod_oc, colour = orig.ident)) + 
    geom_quasirandom(groupOnX = FALSE) + theme_classic() + 
    scale_color_manual(values = c("#5CB85C", "#F0AD4E", "#D9534F")) +
    xlab("Pseudotime") + ylab("Cell types") +  labs(colour = " ") 


# In[232]:


cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
options(repr.plot.width=8, repr.plot.height=4)
ggplot(data.pseudo, aes(monocle3_pseudotime, orig.ident, fill = orig.ident)) + geom_boxplot()


# In[233]:


ggplot(data.pseudo, aes(monocle3_pseudotime, orig.ident, fill = bimod_oc)) + geom_boxplot()


# In[234]:


ggplot(data.pseudo, aes(monocle3_pseudotime, bimod_oc, fill = orig.ident)) + geom_boxplot()


# In[235]:


options(repr.plot.width=8, repr.plot.height=6)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters_pca_RNA, monocle3_pseudotime), fill = seurat_clusters_pca_RNA)) + geom_boxplot()


# In[236]:


cds <- estimate_size_factors(cds)


# In[237]:


options(repr.plot.width=4, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime") 


# In[238]:


options(repr.plot.width=3.5, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "bimod_oc") 


# In[452]:


options(repr.plot.width=7, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "orig.ident") 


# In[449]:


genes = c("NKAIN2","GALNTL6","NTNG1")
lineage_cds = cds[rowData(cds)$gene_short_name %in% genes,]


# In[450]:


cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "bimod_oc") 


# ____

# ## Graph-autocorrelation analysis -- Find genes that change as a function of pseudotime
# 

# In[361]:


cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)


# In[364]:


cds_graph_test_results %>% arrange(q_value) %>% filter(status == "OK") %>% head() 


# In[378]:


Track_genes <- cds_graph_test_results[,c(1,7,4,6,2,1,3)] %>% 
  dplyr::arrange(desc(morans_I),p_value)
head(Track_genes)


# In[379]:


## morans_I arranged

Track_genes <- cds_graph_test_results[,c(1,7,4,6,2,1,3)] %>% 
  dplyr::arrange(desc(morans_I),p_value)
head(Track_genes)


# In[380]:


gene <- Track_genes %>% arrange(q_value) %>% filter(status == "OK") %>% head()


# In[381]:


options(repr.plot.width=8, repr.plot.height=9)
DefaultAssay(sec) <- "RNA"
FeaturePlot(sec, reduction = "umap", 
            features = rownames(gene))


# In[382]:


options(repr.plot.width=20, repr.plot.height=8)
RidgePlot(sec, features = rownames(gene), sort = T)


# In[392]:


options(repr.plot.width=8, repr.plot.height=5)
plot_cells(cds, genes=Track_genes$gene_short_name[1:4],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,)


# ## Find_gene_modules()
# 

# In[407]:


pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))


# In[408]:


Track_genes_2 <- pr_graph_test_res[,c(1,7,4,6,2,1,3)] %>% 
  dplyr::arrange(desc(morans_I),p_value)


# In[409]:


head(Track_genes_2)


# In[410]:


gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))


# In[411]:


table(gene_module_df$module)


# In[412]:


gene_module_df


# In[419]:


saveRDS(gene_module_df,'~/Desktop/AMS_analysis/outputs/s4_monocle_gene_module_df.rds')


# In[413]:


options(repr.plot.width=4, repr.plot.height=10)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=partitions(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


# In[414]:


options(repr.plot.width=8, repr.plot.height=4)

plot_cells(cds, 
           genes=gene_module_df %>% filter(module %in% c(2,9,14,17,4)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)

genes.M=gene_module_df %>% filter(module %in% c(22))
genes.M
write_csv(genes.M, "~/Desktop/XMAS_analysis/monocle3_d56_module9.csv")
# In[420]:


cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=cds$bimod_oc)


# In[422]:


options(repr.plot.width=4, repr.plot.height=10)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=cds$bimod_oc)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


# In[424]:


genes.M <- gene_module_df %>% filter(module %in% c(4))
genes.M
write.csv(genes.M, "~/Desktop/AMS_analysis/outputs/s4_monocle3_FP_module4.csv")


# ## Development related modules

# In[427]:


# 可以将这些在轨迹上变化的pr_deg_ids继续分为小模块
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)


# In[430]:


options(repr.plot.width=5, repr.plot.height=15)

cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$bimod_oc)
agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                    scale="column", clustering_method="ward.D2")


# In[432]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(66,13,15,90,3,57)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


# In[440]:


genes.M <- gene_module_df %>% filter(module %in% c(3,57,90))


# In[441]:


genes.M
write.csv(genes.M, "~/Desktop/AMS_analysis/outputs/s4_monocle3_NbM_modules.csv")


# In[438]:


genes.M <- gene_module_df %>% filter(module %in% c(57))


# In[439]:


genes.M
write.csv(genes.M, "~/Desktop/AMS_analysis/outputs/s4_monocle3_NbM_module57.csv")


# In[444]:


genes = c("NKAIN2","GALNTL6","NTNG1")
lineage_cds = cds[rowData(cds)$gene_short_name %in% genes,]


# In[448]:


options(repr.plot.width=5, repr.plot.height=9)
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="bimod_oc",min_expr=0) + 
    scale_color_brewer(palette = "Paired") +  labs(colour = " ") 


# In[586]:


options(repr.plot.width=5, repr.plot.height=9)
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="orig.ident",min_expr=0) +  labs(colour = " ") 


# In[ ]:





# In[469]:


ml_oc_des <- c(
    "0_Rgl1", 
    "1_ProgFP.G2M",
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
    "20_ProgFP.G2M",
    "21_Rgl3",
    "22_NbM",
    "23_NProg"
)


# In[470]:


ml_oc_des <- gsub(".*_", "", ml_oc_des)


# In[471]:


unique(ml_oc_des)


# In[472]:


cellType_merged_colors <-  viridis(11)
names(cellType_merged_colors) <- c("Rgl1",  "Rgl2", "Rgl3",
                                   "ProgFP.S", "ProgFP.G2M", "ProgFP.Development", "NProg",
                                   "NbM","GabaNb", 
                                   "DA.Neuron projection","DA.Synaptic assembly")


# In[473]:


sec$bimod_oc_des <- ml_oc_des[sec@meta.data$seurat_clusters_pca_RNA_0.3]
sec$bimod_oc_des <- factor(x = sec$bimod_oc_des, levels = names(cellType_merged_colors))


# In[532]:


cds <- subset(sec, subset = bimod_oc %in% c("DA","NbM","NProg","ProgFP") 
              & !(sec$bimod_oc_des %in% c('ProgFP.G2M', 'ProgFP.S','DA.Neuron projection'))
              & !(sec$bimod_oc %in% c("ProgFP","NProg","NbM") & sec$orig.ident == 'New-D28')
              & !(sec$bimod_oc %in% c("ProgFP","NProg") & sec$orig.ident == 'New-D16'))


# In[533]:


cds


# In[534]:


options(repr.plot.width=18, repr.plot.height=4)
p1<-DimPlot(cds, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1,label=TRUE)
p2<-DimPlot(cds, group.by="bimod_oc_des", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1,label=FALSE)
p3<-DimPlot(cds, group.by="bimod_oc", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=TRUE)
p1+p2+p3


# In[535]:


cds <- as.cell_data_set(cds)


# In[536]:


cds <- preprocess_cds(cds)


# In[537]:


cds <- reduce_dimension(cds)


# In[538]:


options(repr.plot.width=5, repr.plot.height=3)
monocle3::plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
monocle3::plot_cells(cds, color_cells_by="bimod_oc", label_cell_groups=FALSE)


# In[539]:


cds <- cluster_cells(cds, resolution = 1e-3)


# In[540]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds,  color_cells_by = "bimod_oc_des", label_groups_by_cluster=FALSE)


# In[541]:


options(repr.plot.width=5, repr.plot.height=3)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster",
           label_cell_groups = TRUE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")


# In[542]:


rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name


# In[543]:


genes <- c("FOXA2","LMX1A","TH","NR4A2","EN1","SOX6")


# In[544]:


# patterns of expression
options(repr.plot.width=10, repr.plot.height=5)
plot_cells(cds,
           genes=genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


# In[545]:


cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)


# In[546]:


get_earliest_principal_node <- function(cds, time_bin="New-D14"){
  cell_ids <- which(colData(cds)[, "orig.ident"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


# In[547]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


# In[548]:


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) %in% c(6)]))


# In[510]:


cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "orig.ident") 


# In[511]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


# In[347]:


options(repr.plot.width=5, repr.plot.height=4)

plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "bimod_oc_des",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")


# In[348]:


pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)


# In[349]:


library(ggbeeswarm)


# In[350]:


library(RColorBrewer)
par(mar=c(3,4,2,2))
options(repr.plot.width=8, repr.plot.height=10)
display.brewer.all()


# In[351]:


options(repr.plot.width=8, repr.plot.height=3)
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = orig.ident, colour = bimod_oc)) + 
    geom_quasirandom(groupOnX = FALSE) + theme_classic() + 
    scale_color_brewer(palette = "Paired")  +
    xlab("Pseudotime") + ylab("Cell types") +  labs(colour = " ") 


# In[352]:


options(repr.plot.width=7, repr.plot.height=3)
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = bimod_oc, colour = orig.ident)) + 
    geom_quasirandom(groupOnX = FALSE) + theme_classic() + 
    scale_color_manual(values = c("#5CB85C", "#F0AD4E", "#D9534F")) +
    xlab("Pseudotime") + ylab("Cell types") +  labs(colour = " ") 


# In[353]:


head(pseudotime(cds), 10)


# In[354]:


cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
options(repr.plot.width=8, repr.plot.height=4)
ggplot(data.pseudo, aes(monocle3_pseudotime, orig.ident, fill = orig.ident)) + geom_boxplot()


# In[355]:


ggplot(data.pseudo, aes(monocle3_pseudotime, orig.ident, fill = bimod_oc)) + geom_boxplot()


# In[356]:


ggplot(data.pseudo, aes(monocle3_pseudotime, bimod_oc, fill = orig.ident)) + geom_boxplot()


# In[357]:


options(repr.plot.width=8, repr.plot.height=6)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters_pca_RNA, monocle3_pseudotime), fill = seurat_clusters_pca_RNA)) + geom_boxplot()


# In[358]:


cds <- estimate_size_factors(cds)


# In[359]:


options(repr.plot.width=4, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime") 


# In[360]:


options(repr.plot.width=3.5, repr.plot.height=10)
cds_subset <- cds[genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "bimod_oc") 

