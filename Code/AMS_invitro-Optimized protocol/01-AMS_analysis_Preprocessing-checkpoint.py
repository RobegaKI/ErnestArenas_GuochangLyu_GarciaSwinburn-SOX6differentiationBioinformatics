#!/usr/bin/env python
# coding: utf-8

# # AMS_analysis_Preprocessing

# ## 1. Loading environments

# ### 1.1 Loading R packages
local({
    r <- getOption("repos")
    r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"
    options(repos=r)
})BiocManager::install('Signac')
# In[3]:


install.packages('~/Downloads/Signac_1.12.0.tar.gz', repos=NULL, type="source")


# In[5]:


library(Seurat)
#library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(dplyr)
library(stringr)
library(harmony)
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
library(pheatmap)
library(grid)
library(biovizBase)
library(ggrepel)
library(viridis)
library(reshape2)
library(hues)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(circlize)
library(scales)
library(hrbrthemes)
library(patchwork)
library(cluster)
library(corrplot)
library(qgraph)
library(plyr)

library(DoubletFinder)
# library(scDblFinder)

# library(girdExtra)


# In[6]:


sessionInfo()


# ### 1.2 Functions definition 
# 

# In[7]:


PrctCellExpringGene <- function(object, genes, group.by = "all"){
    if(group.by == "all"){
        prct = unlist(lapply(genes,calc_helper, object=object))
        result = data.frame(Markers = genes, Cell_proportion = prct)
        return(result)
    }

    else{        
        list = SplitObject(object, group.by)
        factors = names(list)

        results = lapply(list, PrctCellExpringGene, genes=genes)
        for(i in 1:length(factors)){
        results[[i]]$Feature = factors[i]
        }
        combined = do.call("rbind", results)
        return(combined)
    }
}

calc_helper <- function(object,genes){
    counts = object[['RNA']]@counts
    ncells = ncol(counts)
    if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
    }else{return(NA)}
}

NbCellExpringGene <- function(object, genes, group.by = "all"){
    if(group.by == "all"){
        nb = unlist(lapply(genes,calc_helper_nb, object=object))
        result = data.frame(Markers = genes, Cell_nb = nb)
        return(result)
    }

    else{        
        list = SplitObject(object, group.by)
        factors = names(list)

        results = lapply(list, NbCellExpringGene, genes=genes)
        for(i in 1:length(factors)){
        results[[i]]$Feature = factors[i]
        }
        combined = do.call("rbind", results)
        return(combined)
    }
}

calc_helper_nb <- function(object,genes){
    counts = object[['RNA']]@counts
    if(genes %in% row.names(counts)){
    length(counts[genes,][counts[genes,] > 0])
    }else{return(NA)}
}

DoMultiBarHeatmap <- function (object, features = NULL, cells = NULL, group.by = "ident", additional.group.by = NULL, additional.group.sort.by = NULL,  cols.use = NULL, group.bar = TRUE, disp.min = -2.5, disp.max = NULL, slot = "scale.data", assay = NULL, label = TRUE, size = 5.5, hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02, combine = TRUE){

  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
      " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot, " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }

  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ",
        paste(bad.sorts, collapse = ", "))
      }
    }
  }

  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object, save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is.null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]
      if (!is.null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }

    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]

    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }

    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * lines.width), FUN = function(x) {
        return(Seurat:::RandomName(length = 20))
      })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }

      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells

      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells

      group.use <- rbind(group.use, placeholder.groups)

      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }

      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), ncol = ncol(x = data.group), dimnames = list(placeholder.cells, colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }

    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])

    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, disp.min = disp.min, disp.max = disp.max, feature.order = features, cell.order = rownames(x = group.use), group.by = group.use[[i]])

    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }

        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))

        #Overwrite if better value is provided
        if (!is.null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is.null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }

        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])

        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

        plot <- suppressMessages(plot +
          annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
          annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
          coord_cartesian(ylim = c(0, y.max), clip = "off")
        )


        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs

          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]], FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
          plot <- plot + geom_text(stat = "identity", data = label.x.pos, aes_string(label = "group", x = "label.x.pos"), y = y.max + y.max *  0.03 * 0.5, angle = angle, hjust = hjust, size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}


#####################################################
# Function for getting the marker peaks from the Seurat
# object, based on the current ident
#####################################################
getMarkerPeaks <- function(obj, doublets, n_peaks = 100, min_cells = 200){

  n_clusters = length(unique(as.vector(Idents(obj))))
  
  cells <- Cells(obj)
  
  # remove the doublet cells for better marker peaks
  obj <- subset(obj, cells = cells[!cells %in% doublets])
  
  obj@assays$ATAC@data@x <- obj@assays$ATAC@counts@x
  
  # default preporcessing o Signac applied
  obj <- RunTFIDF(obj, verbose = F)
  obj <- FindTopFeatures(obj, min.cutoff = NULL, verbose = F)
  obj <- RunSVD(
    object = obj,
    reduction.key = 'SVD_',
    reduction.name = 'svd',
    seed.use = 1,
    verbose = F
  )
  
  # get the marker peaks between all clusters
  da_peaks <- FindAllMarkers(
    object = obj,
    slot = "data",
    min.pct = 0.1,
    only.pos = TRUE,
    test.use = 'LR',
    min.cells.group = min_cells,
    verbose = F
  )
  
  # only keep the significant peaks
  metadata = lapply(unique(da_peaks$cluster), function(x){
    return(data.frame(da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$cluster == x,]))
  })
  
  names(metadata) <- paste0("Cluster_", unique(da_peaks$cluster))
  
  rm(da_peaks)
  
  meta_peaks = data.frame()
  
  for (i in names(metadata)) {
    
    cur_meta = data.frame(metadata[[i]])
    cur_meta$peaks = row.names(metadata[[i]])
    cur_meta$cluster = i
    row.names(cur_meta) = c()
    
    meta_peaks = rbind(meta_peaks, cur_meta)
    
    rm(cur_meta)
    
  }
  
  # sort the peaks based on the fold changes and signifigance, choose ones with high FC first
  meta_peaks = meta_peaks[order(meta_peaks$avg_logFC, -meta_peaks$p_val_adj, decreasing = T),]
  
  marker_peaks_set = data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("gene", "cluster"))))
  
  # extract 100 marker peaks for each cluster
  while(nrow(marker_peaks_set) < n_peaks*n_clusters) {
    
    temp = meta_peaks[1, c("gene", "cluster")]
    
    marker_peaks_set = rbind(marker_peaks_set, temp)
    
    meta_peaks = meta_peaks[which(meta_peaks$gene != temp$gene),]
    
    if(length(which(marker_peaks_set$cluster == temp$cluster)) == n_peaks){
      meta_peaks = meta_peaks[which(meta_peaks$cluster != temp$cluster),]
    }
    
    rm(temp)
  }
  
  marker_peaks_set = marker_peaks_set[!is.na(marker_peaks_set$gene),]
  
  marker_peaks_set$cluster= factor(marker_peaks_set$cluster, levels = str_sort(unique(marker_peaks_set$cluster), numeric = T))
  marker_peaks_set = marker_peaks_set[order(marker_peaks_set$cluster),]
  
  rm(meta_peaks)
  
  colnames(marker_peaks_set)[1]<- "peaks"
  
  return(marker_peaks_set)
}

############################################
# get the readcount distribution profile for
# the singlet cells of every cluster
############################################
annotateDoublets <- function(obj, marker_peaks, doublets, k = 15){

  # get the read count distributions on the cluster specific marker peaks for all of the cells
  cell.values <- getCellValues(obj, cells = Cells(obj), marker_peaks_set = marker_peaks, doublets = doublets, k = k)
  
  #get the profile for the singlet cells of each cluster
  singlet.profile <- getProfiles(cell.values)
 
  # just annotate the doublets
  doublets <- cell.values %>% 
    subset(doublet == "doublet")
  
  clusters <- doublets %>% dplyr::select(starts_with("Cluster_")) %>% colnames()
  
  # look at every cluster in the data
  a.class <- lapply(clusters, function(clus){
    
    # just select the doublets that have the current cluster
    # as its most dominant cluster
    t.doublets <- doublets %>%
      subset(a.homotypic == clus, select = clusters)
    
    # get the profile of the current cluster
    t.profile <- singlet.profile[,clus]
    
    # calculate the distance of doublet cell distributions to
    # the profile of the current cluster
    t.dist <- apply(t.doublets, 1, function(cell){
      return(dist(rbind(cell, t.profile)))
    })
    
    # fit a GMM using the mclust package with no class limitations
    fit <- Mclust(t.dist, verbose = F)
    t.class <- fit$classification
    
    t.doublets$dist <- t.dist[rownames(t.doublets)]
    t.doublets$class <- t.class[rownames(t.doublets)]
    t.doublets$max.class <- names(which.max(fit[["parameters"]][["mean"]]))
    t.doublets$num.class <- length(unique(fit$classification))
    return(t.doublets[,c("dist", "class", "max.class", "num.class")])
  }) %>% bind_rows()
  
  doublets$class <- a.class[rownames(doublets),"class"]
  doublets$max.class <- a.class[rownames(doublets),"max.class"]
  doublets$num.class <- a.class[rownames(doublets),"num.class"]
  
  doublets <- lapply(rownames(doublets), function(cell){
    probs <- doublets[cell, clusters]
    
    # only classify the doublets that belong to the final class 
    # doublets[cell,"type"] <- ifelse(doublets[cell,"class"]== doublets[cell,"max.class"], "heterotypic", "homotypic")
    doublets[cell,"d.type"] <- ifelse(doublets[cell,"class"]== doublets[cell,"max.class"] & doublets[cell,"num.class"] > 1,
                                    "heterotypic", "homotypic")
    
    # report the top clusters based on type of the doublet; if homotypic just one, if heterotypic top 2 clusters
    doublets[cell, "d.annotation"] <- ifelse(doublets[cell,"d.type"] == "homotypic",
                                             names(probs)[order(probs, decreasing = T)[1]],
                                             paste(sort(names(probs)[order(probs, decreasing = T)[1:2]]), collapse = ".")) 
    
    return(doublets[cell,!(colnames(doublets) %in% c("a.heterotypic", "a.homotypic", "ident", "class", "max.class","num.class")), drop = FALSE])
  }) %>% bind_rows()

  return(doublets)
}

#####################################################
# Function to get the ditsribution of read counts on
# the marker peaks of clusters
#####################################################
getReadCountDistributions <- function(marker_peaks, read_counts){
  
  t_read_counts = data.frame(merge(marker_peaks, 100*read_counts, by.x = "peaks", by.y = 0, sort = F))
  names(t_read_counts) <- sub("^X", "", names(t_read_counts))
  t_read_counts$ids = unclass(t_read_counts$cluster)
  
  probs = matrix(nrow = ncol(read_counts), ncol = length(unique(t_read_counts$cluster)))
  
  row.names(probs) <- colnames(read_counts)
  colnames(probs) <- unique(t_read_counts$cluster)
  
  for (i in 1:nrow(probs)) {
    density_data = vector()
    
    for(j in 1:nrow(t_read_counts))
      density_data = c(density_data, rep(t_read_counts$ids[j], t_read_counts[j,row.names(probs)[i]]))
    
    for(j in 1:ncol(probs)){
      if(j==1)
        probs[i,j] = ecdf(density_data)(j)
      else
        probs[i,j] = (ecdf(density_data)(j) - ecdf(density_data)((j-1)))
    }
    
  }
  
  return(data.frame(probs))
}

#####################################################
# Function to plot the read count distributions of
# the given cells
#####################################################
plotReadCountDistributions <- function(probs, folder_path){
  
  if(!dir.exists(folder_path))
    dir.create(folder_path)
  
  sapply(rownames(probs), function(x){long_probs = gather(as.data.frame(bind_rows(probs[x,])), cell_type, probability, factor_key = TRUE)
  ggplot(data=long_probs, aes(x=cell_type, y=probability)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(x) + xlab("Cluster") + ylab("Score") + ylim(0,1) +
    ggsave(filename = paste0(folder_path, "/Probs_", x, ".pdf"))})
  
}

#####################################################
# Function to add read count distributions for the 
# cells in the provided in te Seurat object with the
# provided marker peaks
#####################################################
getCellValues <- function(obj, cells, marker_peaks_set, doublets, k = 15){
  
  obj@assays$ATAC@data@x <- obj@assays$ATAC@counts@x
  
  obj <- RunTFIDF(obj, verbose = F)
  obj <- FindTopFeatures(obj, min.cutoff = NULL, verbose = F)
  obj <- RunSVD(
    object = obj,
    reduction.key = 'SVD_',
    reduction.name = 'svd',
    seed.use = 1,
    verbose = F
  )
  
  # get the neighbor graph of the cells for aggregation
  obj <- FindNeighbors(obj, reduction = "svd", dims = 1:50, k.param = k+1, verbose = F)
  
  cell_annotations <- lapply(cells, function(cell){
    
    # extract the k neighbors of the current cell
    neighbors <- names(obj@graphs$ATAC_nn[cell, obj@graphs$ATAC_nn[cell,] > 0])
    
    # extract the reads for the cell and its k-nearest neighbors
    reads <- Matrix::as.matrix(subset(obj, cells = neighbors, features = marker_peaks_set$peaks)@assays[["ATAC"]]@counts)
    
    no_clusters <- length(unique(marker_peaks_set$cluster))
    
    results = data.frame(matrix(nrow = 0, ncol = no_clusters+1)) %>%
      `colnames<-`(value = c("cell_id", as.character(unique(marker_peaks_set$cluster))))
    
    results[cell,"cell_id"] = cell
    
    # aggregate the reads for the cell by taking the mean with k-nn
    reads <- data.frame(apply(reads, 1, mean)) %>%
      `colnames<-`(value = cell)
    
    if(colSums(reads) == 0){
      results[,-1] <- 0
      
      results[cell, "a.heterotypic"] <- NA
      results[cell, "a.homotypic"] <- NA
        
      return(results)
    }
    
    # calculate the read count distribution of the cell on the marker peaks
    doublet_probs <- reads %>%
      getReadCountDistributions(marker_peaks_set,.) %>% data.frame()
    
    results[cell, colnames(doublet_probs)] <- doublet_probs
    
    # report the 2 clusters that have the highest score
    results[cell, "a.heterotypic"] <- paste(names(doublet_probs)[order(doublet_probs, decreasing = T)[1:2]], collapse = ".")
    
    #report the cluster with highest score
    results[cell, "a.homotypic"] <- names(which.max(doublet_probs))
    
    return(results)
  }) %>% do.call(rbind, .)
  
  # append the type of the cells
  cell_annotations <- mutate(cell_annotations, doublet = ifelse(cell_id %in% doublets, "doublet", "singlet"))
  row.names(cell_annotations) <- cell_annotations$cell_id
  
  cell_annotations[, "ident"] <- Idents(obj)[rownames(cell_annotations)]
  
  return(cell_annotations)
}

#####################################################
# get the read count distribution profile for the
# singlet cells of every cluster
#####################################################
getProfiles <- function(cell.values){

  # get the singlet cells
  singlets <- cell.values %>%
    subset(doublet == "singlet")
  
  clusters <- singlets %>% dplyr::select(starts_with("Cluster_")) %>% colnames()
  
  # for each cluster in the data find the average profile from singlets
  t.profile <- sapply(clusters, function(clus){
    t.singlets <- singlets %>%
      subset(a.homotypic == clus & paste0("Cluster_", ident) == clus, select = clusters)
    
    t.profile <- apply(t.singlets,2,mean)
    return(t.profile)
  })
  
  return(t.profile)
  
}


# ### 1.3 Set up R parameters 
# 

# In[8]:


#Set up global parameters
OS_path <- "/Users/gclyu07/Desktop/AMS_analysis/"
OS_path_datasets <- paste0(OS_path, "dataset/")
OS_path_inputs  <- paste0(OS_path, "inputs/")
OS_path_outputs <- paste0(OS_path, "outputs/")

seed  <- 1121
options(repr.plot.width=16, repr.plot.height=12)
options(future.globals.maxSize = 8000 * 1024^2)

options(repr.matrix.max.rows=100, repr.matrix.max.cols=100)


# In[9]:


orig.ident_colors <- c("#5CB85C", "#F0AD4E", "#D9534F")
names(orig.ident_colors)  <- c("New-D14", "New-D16", "New-D28")
show_col(orig.ident_colors)


# ## 2. Loading datasets

# ### 2.1 Load Cellranger aggregation

# In[10]:


data_AMS  <- Read10X("~/Desktop/AMS_analysis/dataset/AMS_cellranger_aggr_1/outs/count/filtered_feature_bc_matrix")


# In[11]:


AMS <- CreateSeuratObject(data_AMS)


# In[12]:


AMS


# In[13]:


head(colnames(AMS))


# ### 2.2 Add origin of cell identity

# In[14]:


name_AMS <- seq(1:3)
names(name_AMS) <- c("New-D14", "New-D16", "New-D28")


# In[15]:


AMS$orig.ident <- unlist(lapply(strsplit(rownames(AMS@meta.data),"-"), 
                                 function(x) names(name_AMS[name_AMS == x[2]])))


# In[16]:


AMS$orig.ident <- factor(x = AMS$orig.ident, levels = names(orig.ident_colors))


# In[17]:


table(AMS$orig.ident)


# ## 3. Remove bad cells based on QC metrics

# ### 3.1 QC calculation

# In[18]:


#Add percentage of mt in each cell
AMS$percent.mt <- PercentageFeatureSet(AMS, pattern = "^MT-")


# In[19]:


AMS$log10GenesUMI <- log10(AMS$nFeature_RNA) / log10(AMS$nCount_RNA)


# In[20]:


AMS$mitoRatio <- AMS@meta.data$percent.mt / 100


# ### 3.2 QC visualization

# In[21]:


options(repr.plot.width=30, repr.plot.height=5)
for (sample in 1:length(names(name_AMS))){
    tmp <- subset(AMS, subset = orig.ident == names(name_AMS)[sample])
    Idents(object = tmp) <- "orig.ident"
    print(VlnPlot(object = tmp, features = 
                  c("nCount_RNA",  "percent.mt"), 
                  group.by = "orig.ident", cols=orig.ident_colors, pt.size = 0, ncol=6))
}


# In[22]:


options(repr.plot.width=10, repr.plot.height=5)
RidgePlot(object = AMS, features = 'nCount_RNA', group.by="orig.ident", slot="counts", cols=orig.ident_colors) + xlim(c(0, 20000)) & NoLegend()
RidgePlot(object = AMS, features = 'nFeature_RNA', group.by="orig.ident", slot="counts", cols=orig.ident_colors) + xlim(c(0, 5000)) & NoLegend()


# In[23]:


options(repr.plot.width=15, repr.plot.height=8)
FeatureScatter(AMS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by ="orig.ident", cols=orig.ident_colors)


# In[24]:


options(repr.plot.width=15, repr.plot.height=8)
AMS@meta.data %>% 
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    scale_fill_manual(values=orig.ident_colors) +
    ggtitle("Cell counts")


# In[25]:


AMS@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_density(color="black", alpha = 0.5) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 1000) +
    scale_fill_manual(values=orig.ident_colors) 


# In[26]:


options(repr.plot.width=20, repr.plot.height=15)
AMS@meta.data %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = "glm") +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 600) +
    geom_hline(yintercept = 250) +
    facet_wrap(~orig.ident)


# In[27]:


options(repr.plot.width=15, repr.plot.height=8)
AMS@meta.data %>%
    ggplot(aes(x=log10GenesUMI, color = orig.ident, fill=orig.ident)) +
    geom_density(color='black', alpha = 0.5) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    scale_fill_manual(values=orig.ident_colors)


# ### 3.3 QC cutoffs

# In[28]:


QC_nCount_RNA_max <- c(25000,25000,25000)
QC_nCount_RNA_min <- c(1000,1000,1000)
QC_nFeature_RNA_min <- c(250,250,250)
QC_percentMT_max <- c(20,20,20)


# In[29]:


table(AMS$orig.ident)


# In[30]:


# filter out low quality cells
good_cells <- c()
for (sample in 1:length(names(name_AMS))){
    tmp <- subset(AMS, subset = orig.ident == names(name_AMS)[sample])
    Idents(object = tmp) <- "orig.ident"
    tmp <- subset(
      x = tmp,
      subset = nCount_RNA < QC_nCount_RNA_max[sample] &
        nCount_RNA > QC_nCount_RNA_min[sample] &
        nFeature_RNA > QC_nFeature_RNA_min[sample] &
        percent.mt < QC_percentMT_max[sample]
    )
    good_cells <- c(good_cells, colnames(tmp))
}


# In[31]:


AMS_tmp <- subset(AMS, cells = good_cells)


# In[32]:


table(AMS_tmp$orig.ident)


# In[33]:


options(repr.plot.width=30, repr.plot.height=5)
for (sample in 1:length(names(name_AMS))){
    tmp <- subset(AMS_tmp, subset = orig.ident == names(name_AMS)[sample])
    Idents(object = tmp) <- "orig.ident"
    print(VlnPlot(object = tmp, features = 
                  c("nCount_RNA",  "percent.mt"), 
                  group.by = "orig.ident", cols=orig.ident_colors, ncol = 4, pt.size = 0))
}


# In[34]:


AMS <- AMS_tmp


# ## 4. Doublets detection
# 

# In[35]:


tmp_doublets.list <- SplitObject(AMS, split.by = "orig.ident")


# In[36]:


tmp_doublets.list


# In[37]:


tmp_doublets.list <- lapply(X = tmp_doublets.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
    x <- ScaleData(x, features = rownames(x))
    x <- RunPCA(x)
    x <- FindNeighbors(x, reduction = "pca", dims = 1:15)
    x <- FindClusters(x, resolution = 1)
    x <- RunUMAP(x, reduction = "pca", dims = 1:15)
})

#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[1]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[1]))
print(bcmvn)#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[2]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[2]))
print(bcmvn)#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[3]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[3]))
print(bcmvn)
# In[38]:


best_pk.v <- c(0.005, 0.005, 0.005)
names(best_pk.v) <- names(tmp_doublets.list)


# In[39]:


for (sample in 1:length(tmp_doublets.list)) {

    homotypic.prop <- modelHomotypic(tmp_doublets.list[[sample]]@meta.data$seurat_clusters)
    nExp_poi <- round(0.15*nrow(tmp_doublets.list[[sample]]@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    pN_value <- 0.25

    tmp_doublets.list[[sample]] <- doubletFinder_v3(tmp_doublets.list[[sample]], PCs = 1:15, pN = pN_value, pK = best_pk.v[sample], nExp = nExp_poi.adj, reuse.pANN = FALSE, sct=FALSE)
}


# In[40]:


DF_df <- c()
for (sample in 1:length(tmp_doublets.list)) {
    DF_tmp <- tmp_doublets.list[[sample]]@meta.data[,(ncol(tmp_doublets.list[[sample]]@meta.data)-1):ncol(tmp_doublets.list[[sample]]@meta.data)]
    colnames(DF_tmp) <- c("pANN", "DF.classifications")
    DF_df <- rbind(DF_df, DF_tmp)
}


# In[41]:


AMS@meta.data[c("pANN", "DF.classifications")] <- DF_df[colnames(AMS[["RNA"]]), c("pANN", "DF.classifications")]


# In[42]:


tmp_doublets.list <- NULL


# In[43]:


table(AMS$orig.ident,AMS$DF.classifications)


# In[44]:


AMS <- subset(AMS, DF.classifications == 'Singlet')


# In[45]:


table(AMS$orig.ident,AMS$DF.classifications)


# In[ ]:





# ## 5. Processing RNA assay
# 

# ### 5.1 Normalization, merging & scaling of RNA assay

# In[46]:


AMS.list <- SplitObject(AMS, split.by = "orig.ident")


# In[47]:


AMS.list <- lapply(X = AMS.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# In[48]:


AMS <- merge(AMS.list[[1]], y = AMS.list[2:3], merge.data = TRUE)


# In[49]:


AMS <- FindVariableFeatures(AMS, selection.method = "vst", nfeatures = 2000)


# In[50]:


options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = AMS, selection.method = "vst", log = TRUE)


# In[52]:


AMS <- ScaleData(AMS, features = rownames(AMS))


# In[53]:


library(harmony)
AMS <- NormalizeData(AMS) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
AMS <- RunHarmony(AMS, group.by.vars = "orig.ident")


# In[54]:


AMS


# ### 5.2 Reducing dimension of RNA assay

# In[55]:


AMS <- RunPCA(AMS)


# ### 5.3 Finding clusters in RNA assay
# 

# In[65]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(AMS, reduction = "harmony", ndims = 50)


# In[66]:


AMS <- FindNeighbors(AMS, reduction = "harmony", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA")


# In[67]:


testing_clusters <- AMS
for (clust_res in seq(0.2, 3, by=0.5)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[68]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[69]:


AMS <- FindClusters(AMS, resolution = 1.0, graph.name="snn_pca_RNA")
AMS$seurat_clusters_pca_RNA <- AMS$seurat_clusters
AMS$seurat_clusters <- NULL


# In[70]:


table(AMS$seurat_clusters_pca_RNA)


# ### 5.4 Cells projection in 2D

# In[71]:


AMS <- RunUMAP(AMS, reduction = "harmony", reduction.name="umap_pca_RNA", nn.name="snn_pca_RNA", dims = 1:40)


# In[72]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(AMS, group.by="seurat_clusters_pca_RNA", reduction="umap_pca_RNA", seed=seed, pt.size=1, label=TRUE)


# In[73]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(AMS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[74]:


options(repr.plot.width=9, repr.plot.height=7)
FeaturePlot(AMS, features = "nCount_RNA", reduction = "umap_pca_RNA", cols=c("lightgrey", "red"),pt.size = 0.75, max.cutoff = 3000)
FeaturePlot(AMS, features = "nFeature_RNA", reduction = "umap_pca_RNA", cols=c("lightgrey", "red"), pt.size = 0.75)


# In[75]:


options(repr.plot.width=20, repr.plot.height=12)

FeaturePlot(AMS, features = "nCount_RNA", split.by = 'orig.ident', 
            reduction = "umap_pca_RNA", cols=c("lightgrey", "red"),pt.size = 0.15, max.cutoff = 3000) + patchwork::plot_layout(ncol = 3, nrow = 2)


# In[76]:


options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, features = "nCount_RNA", split.by = "orig.ident", 
            reduction = "umap_pca_RNA", pt.size = 0.75) + patchwork::plot_layout(ncol = 3, nrow = 2)


# In[77]:


options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, features = "nFeature_RNA", split.by = "orig.ident", 
            reduction = "umap_pca_RNA", pt.size = 0.75) + patchwork::plot_layout(ncol = 3, nrow = 2)


# In[78]:


options(repr.plot.width=25, repr.plot.height=6)
DimPlot(AMS,  group.by="orig.ident", split.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, ncol=3, label=FALSE)


# In[79]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(AMS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.2, label=FALSE)


# ### 5.5 Sample proportion
# 

# In[80]:


t <- table(Cluster=AMS$seurat_clusters_pca_RNA, Batch=AMS@meta.data[["orig.ident"]])
t <- t[,rev(names(orig.ident_colors))]
t_percent <- round(prop.table(t, 2) * 100 ,digits = 2)/ncol(t)


# In[81]:


options(repr.plot.width=20, repr.plot.height=12)
barplot(t(t_percent), xlab = "Cluster", ylab="Percentage of cells", 
        legend = TRUE, ylim = c(0, round_any(as.integer(max(rowSums(t_percent)))+5, 5, f = ceiling)), col = rev(orig.ident_colors), args.legend = list(bty = "n", x = "top", ncol = 3))


# In[82]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$orig.ident, AMS@meta.data$seurat_clusters_pca_RNA))
for (row in 1:nrow(mat)) mat[row,] <- as.integer(mat[row,]/(sum(mat[row,])/min(rowSums(mat))))
circos.clear()
par(cex = 0.8)

chordDiagram(mat, annotationTrack = c("grid", "axis"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))/3), big.gap = 20, small.gap = 2, order = c(colnames(mat), names(orig.ident_colors)), grid.col=rev(orig.ident_colors))
circos.track(track.index = 1, panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0.3, 0.5))}, bg.border = NA)


# ## 6. Cell cycle scoring

# In[83]:


AMS <- CellCycleScoring(object=AMS, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)


# In[84]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(AMS, features = c('S.Score','G2M.Score'), reduction = "umap_pca_RNA", cols=c("lightgrey", "blue"), pt.size=1, min.cutoff = 0.05)


# In[85]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(AMS, features = c('S.Score','G2M.Score'), reduction = "umap_pca_RNA", cols=c("purple", "yellow"), pt.size=1, min.cutoff = 0.05)


# In[86]:


cc_colors <- c("#2A75CB", "#F5563D", "#F5AB00")
names(cc_colors)  <- c("S", "G2M", "G1")
show_col(cc_colors)


# In[87]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(AMS, group.by="Phase", reduction="umap_pca_RNA", shuffle=TRUE, seed=seed, cols=cc_colors, pt.size=1, label=FALSE)


# In[88]:


table(AMS$orig.ident,AMS$Phase)


# In[89]:


RidgePlot(AMS, group.by="orig.ident", features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)


# In[90]:


AMS <- ScaleData(AMS, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(AMS))


# In[91]:


AMS <- RunPCA(AMS, npcs = 50, verbose = FALSE)
AMS <- FindNeighbors(AMS, dims = 1:30)
AMS <- FindClusters(AMS, resolution = 1.0)
AMS$seurat_clusters_pca_RNA <- AMS$seurat_clusters


# In[96]:


AMS <- RunUMAP(AMS, reduction = "harmony", reduction.name="umap", nn.name="snn_pca_RNA", dims = 1:50)


# In[97]:


options(repr.plot.width=9, repr.plot.height=7)
DimPlot(AMS, group.by="orig.ident", reduction="umap", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[ ]:





# In[98]:


AMS <- RunUMAP(AMS, reduction = "harmony", reduction.name="umap", nn.name="snn_pca_RNA",min.dist = 0.4, dims = 1:50)


# In[99]:


options(repr.plot.width=9, repr.plot.height=7)
DimPlot(AMS, group.by="orig.ident", reduction="umap", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)
DimPlot(AMS,reduction = "umap",pt.size=1,group.by = "Phase")


# ## 7. Identify clusters

# ### 7.1 Genes markers expression

# In[100]:


#DA
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap", 
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3"))


# In[103]:


#Markers requested by NW
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap", slot = "data", pt.size = 0.1,
            features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"))


# In[104]:


#DA, slot = "scale.data",
VlnPlot(AMS, slot = "scale.data", features  =c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3"),
    pt.size = 0.2, ncol = 2)


# In[105]:


# mesencephalic
options(repr.plot.width=18, repr.plot.height=10)
FeaturePlot(AMS, reduction = "umap", 
            features = c("OTX1", "OTX2", "LMX1A", "EN1", "PITX2","SIM2"), ncol=3)


# In[106]:


options(repr.plot.width=18, repr.plot.height=10)
FeaturePlot(AMS, reduction = "umap", 
            features = c("BARHL1", "HOTAIRM1", "HOXA2", "HOXB2", "GATA3","GBX2"), ncol=3)


# ### 7.2 HVGs in each RNA clusters
# 

# In[107]:


AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)


# In[108]:


AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[109]:


AMS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7


# In[110]:


top7


# In[111]:


topDiffCluster <- AMS.markers_best %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


# In[112]:


sample_df <- as.data.frame(lapply(split(topDiffCluster, topDiffCluster$cluster), function(x) c(x$gene,rep("None", 100-length(x$gene)))))
colnames(sample_df) <- levels(topDiffCluster$cluster)
sample_df


# In[113]:


#Use RNA assay to look at gene expression across cluster, but integrated assay to look at genes differences driving clustering
#https://github.com/satijalab/seurat/issues/1717
options(repr.plot.width=25, repr.plot.height=25)
DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "seurat_clusters_pca_RNA"
DoMultiBarHeatmap(
    subset(AMS, downsample = 200),
    features = top7$gene,
    cells = NULL,
    group.by = "seurat_clusters_pca_RNA",
    additional.group.by = c("orig.ident"),
    additional.group.sort.by = c("orig.ident"),
    cols.use = list(orig.ident=orig.ident_colors),
    group.bar = TRUE,
    disp.min = -2.5,
    disp.max = NULL,
    slot = "scale.data",
    assay = "RNA",
    label = TRUE,
    size = 5.5,
    hjust = 0,
    angle = 45,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    group.bar.height = 0.02,
    combine = TRUE
    ) + theme(text = element_text(size = 8))


# In[114]:


write.csv(sample_df,"~/Desktop/AMS_analysis/0725_AMS_clstr_Top100HVGs.csv",row.names = FALSE)


# ## 8. Cell type annotation

# In[115]:


mid <- readRDS("~/Desktop/XMAS_analysis/SL_midbrain_2.rds")


# In[116]:


mid <- SCTransform(mid)
DefaultAssay(AMS) <- "RNA"

Am_anchor  <- FindTransferAnchors(reference = mid, query = AMS, 
                                  normalization.method = 'SCT', recompute.residuals = TRUE,
                                  reduction = 'cca', dims = 1:10)
Am_predictions <- TransferData(anchorset = Am_anchor, refdata = mid$LRprediction_labels, dims = 1:10,
                              weight.reduction = 'cca')
AMS <- AddMetaData(AMS, Am_predictions$predicted.id, col.name = 'Cell_type_SL')
for (prediction_score in colnames(Am_predictions)[!colnames(Am_predictions) %in% c("predicted.id", "prediction.score.max")]){
  AMS <- AddMetaData(AMS, Am_predictions[prediction_score], col.name = paste0("Cell_type_SL_",prediction_score))
}


# In[117]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="Cell_type_SL", reduction = "umap_pca_RNA", pt.size=1, label=TRUE) + NoAxes() + ggtitle("RNA")


# In[118]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="orig.ident", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE, cols=orig.ident_colors ) + NoAxes() + ggtitle("Original Identity")
DimPlot(AMS, group.by="Cell_type_SL", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("RNA")


# In[119]:


table(AMS$Cell_type_SL)


# In[120]:


table(AMS$orig.ident, AMS$Cell_type_SL)


# In[121]:


colors <- hue_pal(h.start = 30)(36)

expression <- as.data.frame(as.matrix(AMS@meta.data[,colnames(AMS@meta.data)[colnames(AMS@meta.data) %like% "Cell_type_SL_prediction.score."]]))
colnames(expression) <- gsub("Cell_type_SL_prediction.score.", "", colnames(expression))

expression <- expression[,c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")]
    expression <-cbind(expression,Clusters=AMS$seurat_clusters_pca_RNA)
expression_melt <- reshape2::melt(expression,id=c("Clusters"))
colnames(expression_melt)[3] <- "prediction.score"

expression_melt <- expression_melt %>% mutate(Clusters=factor(Clusters,levels=levels(AMS$seurat_clusters_pca_RNA)))


# In[122]:


options(repr.plot.width=15, repr.plot.height=15)
p <- ggplot(expression_melt, aes(x=Clusters, y=prediction.score,fill=Clusters))
p +geom_violin(scale="width") + geom_boxplot(width=0.1,outlier.shape = NA,position=position_dodge(1),fill="white")+
  theme_bw()+scale_fill_manual(values=colors)+ylim(0,1)+  
theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="right") +facet_grid(variable ~ .,scales="free") +theme(strip.text.y = element_text(angle = 0))


# In[123]:


desired_order <- c("D11", "D16", "D28", "D42","D56")
cell_class <- c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","
OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")

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

percentages <- tapply(AMS$Cell_type_SL, AMS@meta.data[["orig.ident"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

A_labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[124]:


A_labels #percentage > 3%


# In[125]:


options(repr.plot.width=17, repr.plot.height=3)

n_colors <- 15
color_palette <- viridis::viridis(n_colors, option = "D")

desired_order <- c("New-D14", "New-D16", "New-D28")

t <- NULL
ggplot(t, aes(x = factor(AMS@meta.data[["orig.ident"]], levels = desired_order), fill = AMS$Cell_type_SL)) +
  geom_bar(position = "fill", color='white',width=0.9) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(limits = desired_order) +
  scale_y_continuous(labels = scales::percent) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_blank()
  ) 




# In[126]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL' == 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# In[127]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$orig.ident, AMS@meta.data$Cell_type_SL))
mat <- mat[,c("DA", "DA0", "Gaba", "GabaNb", "NbM",  "NbML1", "NProg", "OMTN", "ProgBP", "ProgFP", "Rgl1", "Rgl2", "Rgl3","RN")]
circos.clear()
par(cex = 1)
chordDiagram(mat, big.gap = 20, small.gap = 2, order = c(colnames(mat), rownames(mat)), grid.col = orig.ident_colors)  


# In[128]:


options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(AMS, slot = 'scale.data', group.by = "Cell_type_SL", features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"),ncol=2)


# In[129]:


saveRDS(AMS,"~/Desktop/AMS_analysis/0725_AMS_AFTER_LT.rds")


# ### 8.1 SL developing midbrain dataset (lt,abv 0.3)

# In[131]:


t <- ifelse(Am_predictions$prediction.score.max > 0.3,
        Am_predictions$predicted.id,
            'Unknown')


# In[132]:


AMS <- AddMetaData(AMS, t, col.name = 'Cell_type_SL_abv0.3')


# In[133]:


options(repr.plot.width=12, repr.plot.height=8)
DimPlot(AMS, group.by="Cell_type_SL_abv0.3", reduction = "umap_pca_RNA", pt.size=1, label=TRUE)


# In[134]:


table(AMS$orig.ident, AMS$Cell_type_SL_abv0.3)


# In[135]:


desired_order <- c("D11", "D16", "D28", "D42","D56")
cell_class <- c("Unknown","DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")

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

percentages <- tapply(AMS$Cell_type_SL_abv0.3, AMS@meta.data[["orig.ident"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

A_labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 0, paste0(round(p, 1), "%"), " ")))


# In[136]:


A_labels


# ### 8.2 Subsetting with thresholds as 0.3

# In[137]:


non_un <- subset(AMS, subset = Cell_type_SL_abv0.3 != "Unknown")


# In[138]:


non_un
table(non_un$orig.ident)


# In[140]:


options(repr.plot.width=17, repr.plot.height=3)

n_colors <- 15
color_palette <- viridis::viridis(n_colors, option = "D")
desired_order <- c("New-D14", "New-D16", "New-D28")

t <- NULL

ggplot(t, aes(x = factor(non_un@meta.data[["orig.ident"]], levels = desired_order), fill = non_un$Cell_type_SL_abv0.3)) +
  geom_bar(position = "fill", color='white',width=0.9) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(limits = desired_order) +
  scale_y_continuous(labels = scales::percent) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_blank()
  )


# In[141]:


desired_order <- c("New-D14", "New-D16", "New-D28")
cell_class <- c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")

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

percentages <- tapply(non_un$Cell_type_SL_abv0.3, non_un@meta.data[["orig.ident"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 0, paste0(round(p, 1), "%"), " ")))
                               
labels


# In[142]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL_abv0.3' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL_abv0.3' == 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL_abv0.3' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# ## 9. Re-projecting subsets in RNA

# In[143]:


table(non_un$orig.ident)


# In[146]:


non_un <- RunPCA(non_un)
options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(non_un, reduction = "harmony", ndims = 50)


# In[147]:


non_un <- FindNeighbors(non_un, reduction = "harmony", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA_0.3")


# In[148]:


testing_clusters <- non_un
for (clust_res in seq(0.4, 3, by=0.3)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA_0.3")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[149]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_0.3_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[151]:


non_un <- FindClusters(non_un, resolution = 1.3, graph.name="snn_pca_RNA_0.3")
non_un$seurat_clusters_pca_RNA_0.3 <- non_un$seurat_clusters
non_un$seurat_clusters <- NULL


# In[152]:


non_un <- RunUMAP(non_un, reduction = "harmony", reduction.name="umap_pca_RNA_0.3", nn.name="snn_pca_RNA_0.3", dims = 1:40, min.dist = 0.5)


# In[153]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=TRUE)


# In[154]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="Cell_type_SL_abv0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=FALSE)


# In[155]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="orig.ident", reduction="umap_pca_RNA_0.3", cols=orig.ident_colors, seed=seed, pt.size=1, label=FALSE)


# In[156]:


#DA, slot = "scale.data",
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(non_un, reduction = "umap_pca_RNA_0.3", slot = "data", pt.size = 0.1,
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3"))


# In[157]:


#Markers requested by NW
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(non_un, reduction = "umap_pca_RNA_0.3", slot = "data", pt.size = 0.1,
            features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"))


# In[159]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(non_un, reduction = "umap_pca_RNA_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL_abv0.3' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(non_un, reduction = "umap_pca_RNA_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL_abv0.3' == 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(non_un, reduction = "umap_pca_RNA_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL_abv0.3' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# _______

# In[160]:


saveRDS(non_un,"/Users/gclyu07/Desktop/AMS_analysis/0725_s2_AMS_ABV0.3_FILTERED.rds")


# ## 10. HVGs in each cell type (abv 0.3)

# In[178]:


AMS <- non_un


# In[179]:


DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "Cell_type_SL_abv0.3"


# In[ ]:


AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)


# In[ ]:


AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[ ]:


AMS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7


# In[ ]:


top7


# In[ ]:


topDiffCluster <- AMS.markers_best %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


# In[ ]:


sample_df <- as.data.frame(lapply(split(topDiffCluster, topDiffCluster$cluster), function(x) c(x$gene,rep("None", 100-length(x$gene)))))
colnames(sample_df) <- levels(topDiffCluster$cluster)
sample_df


# In[ ]:


#Use RNA assay to look at gene expression across cluster, but integrated assay to look at genes differences driving clustering
#https://github.com/satijalab/seurat/issues/1717
options(repr.plot.width=25, repr.plot.height=25)
DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "Cell_type_SL"
DoMultiBarHeatmap(
    subset(AMS, downsample = 200),
    features = top7$gene,
    cells = NULL,
    group.by = "Cell_type_SL_abv0.3",
    additional.group.by = c("orig.ident"),
    additional.group.sort.by = c("orig.ident"),
    cols.use = list(orig.ident=orig.ident_colors),
    group.bar = TRUE,
    disp.min = -2.5,
    disp.max = NULL,
    slot = "scale.data",
    assay = "RNA",
    label = TRUE,
    size = 5.5,
    hjust = 0,
    angle = 45,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    group.bar.height = 0.02,
    combine = TRUE
    ) + theme(text = element_text(size = 8))


# In[ ]:


write.csv(sample_df,paste0(OS_path_outputs,"0725_AMS_Top100HVGs_Celltype_abv0.3.csv"),row.names = FALSE)


# ## 11. HVGs in each cluster (abv 0.3)

# In[ ]:


AMS$seurat_clusters_pca_RNA


# In[ ]:


DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "seurat_clusters_pca_RNA_0.3"


# In[ ]:


AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)


# In[ ]:


AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[ ]:


AMS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7


# In[ ]:


top7


# In[ ]:


topDiffCluster <- AMS.markers_best %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


# In[ ]:


sample_df <- as.data.frame(lapply(split(topDiffCluster, topDiffCluster$cluster), function(x) c(x$gene,rep("None", 100-length(x$gene)))))
colnames(sample_df) <- levels(topDiffCluster$cluster)
sample_df


# In[ ]:


#Use RNA assay to look at gene expression across cluster, but integrated assay to look at genes differences driving clustering
#https://github.com/satijalab/seurat/issues/1717
options(repr.plot.width=25, repr.plot.height=25)
DefaultAssay(AMS) <- "RNA"
Idents(object = AMS) <- "seurat_clusters_pca_RNA_0.3"
DoMultiBarHeatmap(
    subset(AMS, downsample = 200),
    features = top7$gene,
    cells = NULL,
    group.by = "seurat_clusters_pca_RNA_0.3",
    additional.group.by = c("orig.ident"),
    additional.group.sort.by = c("orig.ident"),
    cols.use = list(orig.ident=orig.ident_colors),
    group.bar = TRUE,
    disp.min = -2.5,
    disp.max = NULL,
    slot = "scale.data",
    assay = "RNA",
    label = TRUE,
    size = 5.5,
    hjust = 0,
    angle = 45,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    group.bar.height = 0.02,
    combine = TRUE
    ) + theme(text = element_text(size = 8))


# In[ ]:


write.csv(sample_df,paste0(OS_path_outputs,"0725_AMS_Top100HVGs_clstr_abv0.3.csv"),row.names = FALSE)

