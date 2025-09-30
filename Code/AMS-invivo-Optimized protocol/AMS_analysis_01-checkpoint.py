#!/usr/bin/env python
# coding: utf-8

# # AMS_Single Cell Gene Expression 01
# 

# ## 1. Loading environments

# ### 1.1 Loading R packages

# In[2]:


library(Seurat)
library(Signac)
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


# In[2]:


sessionInfo()


# ### 1.2 Functions definition 
# 

# In[3]:


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

# In[4]:


#Set up global parameters
OS_path <- "/Users/gclyu07/Desktop/R_data/"
OS_path_datasets <- paste0(OS_path, "dataset/")
OS_path_inputs  <- paste0(OS_path, "inputs/")
OS_path_outputs <- paste0(OS_path, "outputs/")

seed  <- 1121
options(repr.plot.width=16, repr.plot.height=12)
options(future.globals.maxSize = 8000 * 1024^2)

options(repr.matrix.max.rows=100, repr.matrix.max.cols=100)


# In[136]:


orig.ident_colors <- c("#5CB85C", "#D9534F","#428BCA")
names(orig.ident_colors)  <- c("795011","794999","795008")
show_col(orig.ident_colors)


# ## 2. Loading datasets

# ### 2.1 Load Cellranger aggregation

# In[147]:


data_AMS  <- Read10X(paste0(OS_path_datasets, "cellranger_aggr_wo134/R_cellranger_aggr/outs/count/filtered_feature_bc_matrix"))


# In[225]:


AMS <- CreateSeuratObject(data_AMS)


# In[226]:


AMS


# In[227]:


head(colnames(AMS))


# ### 2.2 Add origin of cell identity

# In[228]:


name_AMS <- seq(1:3)
names(name_AMS) <- c("795011","794999","795008")


# In[229]:


AMS$orig.ident <- unlist(lapply(strsplit(rownames(AMS@meta.data),"-"), 
                                 function(x) names(name_AMS[name_AMS == x[2]])))


# In[230]:


AMS$orig.ident <- factor(x = AMS$orig.ident, levels = names(orig.ident_colors))


# In[231]:


table(AMS$orig.ident)


# ### 2.3 Species Identification

# In[232]:


head(rownames(AMS))


# In[233]:


grch38_genes <- grep("^GRCh38-", rownames(AMS), value = TRUE)
grch38_counts <- AMS@assays$RNA@counts[grch38_genes, ]
grch38_positive_counts <- colSums(grch38_counts)
AMS$Human_transcripts <- grch38_positive_counts


# In[234]:


head(grch38_positive_counts)
head(AMS$Human_transcripts)


# In[235]:


grcm39_genes <- grep("^GRCm39-", rownames(AMS), value = TRUE)
grcm39_counts <- AMS@assays$RNA@counts[grcm39_genes, ]
grcm39_positive_counts <- colSums(grcm39_counts)
AMS$Mouse_transcripts <- grcm39_positive_counts


# In[236]:


head(grcm39_positive_counts)
head(AMS$Human_transcripts)


# In[237]:


AMS$Human_Mouse_Ratio <- AMS$Human_transcripts / AMS$Mouse_transcripts


# In[238]:


table(AMS$Human_Mouse_Ratio>1)


# In[239]:


options(repr.plot.width=8, repr.plot.height=6)
FeatureScatter(AMS, feature1 = "Human_transcripts", feature2 = "Mouse_transcripts", group.by ="orig.ident", cols=orig.ident_colors)


# In[240]:


VlnPlot(AMS, feature = "Human_Mouse_Ratio", group.by ="orig.ident", cols=orig.ident_colors)


# In[241]:


AMS_h <- subset(AMS, subset = Human_Mouse_Ratio >1)


# In[242]:


grch38_genes <- grep("^GRCh38-", rownames(AMS), value = TRUE)
AMS_h <- subset(AMS_h, features = grch38_genes)


# In[243]:


head(rownames(AMS_h))


# In[244]:


grch38_genes_2 <- gsub("^GRCh38-", "", grch38_genes)


# In[245]:


head(grch38_genes_2)


# In[248]:


length(rownames(AMS_h))
length(rownames(AMS_h@assays$RNA@counts))
length(rownames(AMS_h@assays$RNA@data))
length(rownames(AMS_h@assays$RNA@scale.data))
length(grch38_genes_2)


# In[252]:


rownames(AMS_h@assays$RNA@counts) <- grch38_genes_2
rownames(AMS_h@assays$RNA@data) <- grch38_genes_2
rownames(AMS_h@assays$RNA@meta.features) <- grch38_genes_2


# In[254]:


head(rownames(AMS_h@assays$RNA@counts))
head(rownames(AMS_h@assays$RNA@data))
head(rownames(AMS_h@assays$RNA@meta.features))
head(grch38_genes_2)


# In[255]:


AMS_h$Human_features <- colSums(AMS_h@assays$RNA@counts>0)


# In[256]:


AMS <- AMS_h


# ## 3. Remove bad cells based on QC metrics

# ### 3.1 QC calculation

# In[257]:


#Add percentage of mt in each cell
AMS$percent.mt <- PercentageFeatureSet(AMS, pattern = "^MT-")


# In[258]:


AMS$log10GenesUMI <- log10(AMS$nFeature_RNA) / log10(AMS$nCount_RNA)


# In[259]:


AMS$mitoRatio <- AMS@meta.data$percent.mt / 100


# ### 3.2 QC visualization

# In[260]:


library(patchwork)


# In[261]:


options(repr.plot.width=30, repr.plot.height=5)
for (sample in 1:length(names(name_AMS))){
    tmp <- subset(AMS, subset = orig.ident == names(name_AMS)[sample])
    Idents(object = tmp) <- "orig.ident"
    print(VlnPlot(object = tmp, features = 
                  c("Human_transcripts",  "percent.mt"), 
                  group.by = "orig.ident", cols=orig.ident_colors, pt.size = 0, ncol=6))
}


# In[262]:


options(repr.plot.width=10, repr.plot.height=5)
RidgePlot(object = AMS, features = 'Human_transcripts', group.by="orig.ident", slot="counts", cols=orig.ident_colors) + xlim(c(0, 20000)) & NoLegend()
RidgePlot(object = AMS, features = 'Human_features', group.by="orig.ident", slot="counts", cols=orig.ident_colors) + xlim(c(0, 5000)) & NoLegend()


# In[263]:


options(repr.plot.width=15, repr.plot.height=8)
FeatureScatter(AMS, feature1 = "Human_transcripts", feature2 = "Human_features", group.by ="orig.ident", cols=orig.ident_colors)


# In[264]:


options(repr.plot.width=15, repr.plot.height=8)
AMS@meta.data %>% 
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    scale_fill_manual(values=orig.ident_colors) +
    ggtitle("Cell counts")


# In[265]:


AMS@meta.data %>% 
    ggplot(aes(color=orig.ident, x=Human_transcripts, fill= orig.ident)) + 
    geom_density(color="black", alpha = 0.5) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 1000) +
    scale_fill_manual(values=orig.ident_colors) 


# In[266]:


options(repr.plot.width=10, repr.plot.height=6)
AMS@meta.data %>% 
    ggplot(aes(x=Human_transcripts, y=Human_features, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = "glm") +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 600) +
    geom_hline(yintercept = 250) +
    facet_wrap(~orig.ident)


# In[267]:


options(repr.plot.width=15, repr.plot.height=8)
AMS@meta.data %>%
    ggplot(aes(x=log10GenesUMI, color = orig.ident, fill=orig.ident)) +
    geom_density(color='black', alpha = 0.5) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    scale_fill_manual(values=orig.ident_colors)


# ### 3.3 QC cutoffs

# In[268]:


QC_nCount_RNA_max <- c(25000,25000,25000)
QC_nCount_RNA_min <- c(1000,1000,1000)
QC_nFeature_RNA_min <- c(250,250,250)
QC_percentMT_max <- c(10,10,10)


# In[269]:


table(AMS$orig.ident)


# In[270]:


# filter out low quality cells
good_cells <- c()
for (sample in 1:length(names(name_AMS))){
    tmp <- subset(AMS, subset = orig.ident == names(name_AMS)[sample])
    Idents(object = tmp) <- "orig.ident"
    tmp <- subset(
      x = tmp,
      subset = Human_transcripts < QC_nCount_RNA_max[sample] &
        Human_transcripts > QC_nCount_RNA_min[sample] &
        Human_features > QC_nFeature_RNA_min[sample] &
        percent.mt < QC_percentMT_max[sample]
    )
    good_cells <- c(good_cells, colnames(tmp))
}


# In[271]:


AMS_tmp <- subset(AMS, cells = good_cells)


# In[272]:


table(AMS_tmp$orig.ident)


# In[273]:


options(repr.plot.width=30, repr.plot.height=5)
for (sample in 1:length(names(name_AMS))){
    tmp <- subset(AMS_tmp, subset = orig.ident == names(name_AMS)[sample])
    Idents(object = tmp) <- "orig.ident"
    print(VlnPlot(object = tmp, features = 
                  c("nCount_RNA",  "percent.mt"), 
                  group.by = "orig.ident", cols=orig.ident_colors, ncol = 4, pt.size = 0))
}


# In[274]:


AMS <- AMS_tmp


# In[275]:


saveRDS(AMS, paste0(OS_path_outputs, "01_AMS_afterQC.rds"))


# ## 4. Doublets detection
# 

# In[290]:


tmp_doublets.list <- SplitObject(AMS, split.by = "orig.ident")


# In[291]:


tmp_doublets.list


# In[292]:


tmp_doublets.list <- lapply(X = tmp_doublets.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
    x <- ScaleData(x, features = rownames(x))
    x <- RunPCA(x)
    x <- FindNeighbors(x, reduction = "pca", dims = 1:15)
    x <- FindClusters(x, resolution = 1)
    x <- RunUMAP(x, reduction = "pca", dims = 1:15)
})


# In[294]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[1]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[1]))
print(bcmvn)


# In[293]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[2]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[2]))
print(bcmvn)


# In[295]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[3]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[3]))
print(bcmvn)


# In[296]:


best_pk.v <- c(0.03, 0.01, 0.01)
names(best_pk.v) <- names(tmp_doublets.list)


# In[297]:


for (sample in 1:length(tmp_doublets.list)) {

    homotypic.prop <- modelHomotypic(tmp_doublets.list[[sample]]@meta.data$seurat_clusters)
    nExp_poi <- round(0.15*nrow(tmp_doublets.list[[sample]]@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    pN_value <- 0.25

    tmp_doublets.list[[sample]] <- doubletFinder_v3(tmp_doublets.list[[sample]], PCs = 1:15, pN = pN_value, pK = best_pk.v[sample], nExp = nExp_poi.adj, reuse.pANN = FALSE, sct=FALSE)
}


# In[298]:


DF_df <- c()
for (sample in 1:length(tmp_doublets.list)) {
    DF_tmp <- tmp_doublets.list[[sample]]@meta.data[,(ncol(tmp_doublets.list[[sample]]@meta.data)-1):ncol(tmp_doublets.list[[sample]]@meta.data)]
    colnames(DF_tmp) <- c("pANN", "DF.classifications")
    DF_df <- rbind(DF_df, DF_tmp)
}


# In[299]:


AMS@meta.data[c("pANN", "DF.classifications")] <- DF_df[colnames(AMS[["RNA"]]), c("pANN", "DF.classifications")]


# In[300]:


tmp_doublets.list <- NULL


# In[301]:


table(AMS$orig.ident,AMS$DF.classifications)


# In[419]:


table(AMS$orig.ident)


# In[361]:


AMS <- subset(AMS, DF.classifications == 'Singlet')


# ## 5. Processing RNA assay
# 

# ### 5.1 Normalization, merging & scaling of RNA assay

# In[362]:


AMS.list <- SplitObject(AMS, split.by = "orig.ident")


# In[363]:


AMS.list <- lapply(X = AMS.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# In[364]:


AMS <- merge(AMS.list[[1]], y = AMS.list[2:3], merge.data = TRUE)


# In[365]:


AMS <- FindVariableFeatures(AMS, selection.method = "vst", nfeatures = 2000)


# In[366]:


options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = AMS, selection.method = "vst", log = TRUE)


# In[367]:


AMS <- ScaleData(AMS, features = rownames(AMS))


# In[368]:


saveRDS(AMS, paste0(OS_path_outputs, "02_AMS_afterScaling.rds"))


# ### 5.2 Reducing dimension of RNA assay

# In[369]:


AMS <- RunPCA(AMS)


# ### 5.3 Finding clusters in RNA assay
# 

# In[370]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(AMS, reduction = "pca", ndims = 50)


# In[371]:


AMS <- FindNeighbors(AMS, reduction = "pca", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA")


# In[372]:


testing_clusters <- AMS
for (clust_res in seq(0.2, 3, by=0.5)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[373]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[374]:


AMS <- FindClusters(AMS, resolution = 1.0, graph.name="snn_pca_RNA")
AMS$seurat_clusters_pca_RNA <- AMS$seurat_clusters
AMS$seurat_clusters <- NULL


# In[89]:


AMS


# In[375]:


table(AMS$seurat_clusters_pca_RNA)


# ### 5.4 Cells projection in 2D

# In[400]:


AMS <- RunUMAP(AMS, reduction = "pca", reduction.name="umap_pca_RNA", nn.name="snn_pca_RNA", dims = 1:40, min.dist = 0.75)


# In[401]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(AMS, group.by="seurat_clusters_pca_RNA", reduction="umap_pca_RNA", seed=seed, pt.size=1, label=TRUE)


# In[402]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(AMS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[403]:


options(repr.plot.width=9, repr.plot.height=7)
FeaturePlot(AMS, features = "Human_transcripts", reduction = "umap_pca_RNA", cols=c("lightgrey", "red"),pt.size = 0.75, max.cutoff = 3000)
FeaturePlot(AMS, features = "Human_features", reduction = "umap_pca_RNA", cols=c("lightgrey", "red"), pt.size = 0.75)


# In[404]:


options(repr.plot.width=20, repr.plot.height=6)

FeaturePlot(AMS, features = "Human_transcripts", split.by = 'orig.ident', 
            reduction = "umap_pca_RNA", cols=c("lightgrey", "red"),pt.size = 0.15, max.cutoff = 3000) + patchwork::plot_layout(ncol = 3, nrow = 1)


# In[405]:


options(repr.plot.width=20, repr.plot.height=6)
FeaturePlot(AMS, features = "Human_features", split.by = "orig.ident", 
            reduction = "umap_pca_RNA", pt.size = 0.75) + patchwork::plot_layout(ncol = 3, nrow = 1)


# In[406]:


options(repr.plot.width=25, repr.plot.height=6)
DimPlot(AMS,  group.by="orig.ident", split.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, ncol=3, label=FALSE)


# In[407]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(AMS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.2, label=FALSE)


# ### 5.5 Sample proportion
# 

# In[408]:


t <- table(Cluster=AMS$seurat_clusters_pca_RNA, Batch=AMS@meta.data[["orig.ident"]])
t <- t[,rev(names(orig.ident_colors))]
t_percent <- round(prop.table(t, 2) * 100 ,digits = 2)/ncol(t)


# In[409]:


options(repr.plot.width=20, repr.plot.height=12)
barplot(t(t_percent), xlab = "Cluster", ylab="Percentage of cells", 
        legend = TRUE, ylim = c(0, round_any(as.integer(max(rowSums(t_percent)))+5, 5, f = ceiling)), col = rev(orig.ident_colors), args.legend = list(bty = "n", x = "top", ncol = 3))


# In[410]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$orig.ident, AMS@meta.data$seurat_clusters_pca_RNA))
for (row in 1:nrow(mat)) mat[row,] <- as.integer(mat[row,]/(sum(mat[row,])/min(rowSums(mat))))
circos.clear()
par(cex = 0.8)

chordDiagram(mat, annotationTrack = c("grid", "axis"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))/3), big.gap = 20, small.gap = 2, order = c(colnames(mat), names(orig.ident_colors)), grid.col=rev(orig.ident_colors))
circos.track(track.index = 1, panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0.3, 0.5))}, bg.border = NA)


# ## 6. Cell cycle scoring

# In[411]:


AMS <- CellCycleScoring(object=AMS, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)


# In[412]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(AMS, features = c('S.Score','G2M.Score'), reduction = "umap_pca_RNA", cols=c("lightgrey", "blue"), pt.size=1, min.cutoff = 0.05)


# In[413]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(AMS, features = c('S.Score','G2M.Score'), reduction = "umap_pca_RNA", cols=c("purple", "yellow"), pt.size=1, min.cutoff = 0.05)


# In[414]:


cc_colors <- c("#2A75CB", "#F5563D", "#F5AB00")
names(cc_colors)  <- c("S", "G2M", "G1")
show_col(cc_colors)


# In[415]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(AMS, group.by="Phase", reduction="umap_pca_RNA", shuffle=TRUE, seed=seed, cols=cc_colors, pt.size=1, label=FALSE)


# In[416]:


table(AMS$orig.ident,AMS$Phase)


# In[417]:


RidgePlot(AMS, group.by="orig.ident", features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)


# ## 7. Identify clusters

# ### 7.1 Genes markers expression

# In[413]:


options(repr.plot.width=10, repr.plot.height=4)
FeaturePlot(AMS, reduction = "umap_pca_RNA", 
            features = c("MCM6", "MKI67"), pt.size=0.05)


# In[567]:


#VLMC,OPC
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap_pca_RNA", 
            features = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "FBLN1", "PDGFRA", "SOX10", 
                         "OLIG1", "MYRF", "CNP","PLP1"), pt.size=0.05)


# In[92]:


#VLMC,OPC
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot_scCustom(AMS, reduction = "umap_pca_RNA", 
            features = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "FBLN1", "PDGFRA", "SOX10", 
                         "OLIG1", "MYRF", "CNP","PLP1"), pt.size=0.05, colors_use = viridis_dark_high)


# In[ ]:





# In[444]:


#DA
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap_pca_RNA", 
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1"), pt.size=0.05)


# In[418]:


#DA
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap_pca_RNA", 
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3"))


# In[419]:


#DA, slot = "scale.data",
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap_pca_RNA", slot = "scale.data", pt.size = 0.1,
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3"))


# In[420]:


#Markers requested by NW
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(AMS, reduction = "umap_pca_RNA", slot = "scale.data", pt.size = 0.1,
            features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"))


# In[685]:


#DA, slot = "scale.data",
options(repr.plot.width=20, repr.plot.height=25)
VlnPlot(AMS, slot = "scale.data", features  =c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3",'PITX3','ALDH1A1','GAD1','GAD2'),
    pt.size = 0.2, ncol = 2)


# In[422]:


# mesencephalic
options(repr.plot.width=18, repr.plot.height=10)
FeaturePlot(AMS, reduction = "umap_pca_RNA", 
            features = c("OTX1", "OTX2", "LMX1A", "EN1", "PITX2","SIM2"), ncol=3)


# In[423]:


options(repr.plot.width=18, repr.plot.height=10)
FeaturePlot(AMS, reduction = "umap_pca_RNA", 
            features = c("BARHL1", "HOTAIRM1", "HOXA2", "HOXB2", "GATA3","GBX2"), ncol=3)


# ### 7.2 HVGs in each RNA clusters
# 

# In[57]:


AMS.markers <- FindAllMarkers(AMS, slot = "data", min.pct = 0.05, logfc.threshold = 0.5, only.pos = TRUE)


# In[58]:


AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[59]:


AMS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7


# In[60]:


top7


# In[61]:


topDiffCluster <- AMS.markers_best %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


# In[62]:


sample_df <- as.data.frame(lapply(split(topDiffCluster, topDiffCluster$cluster), function(x) c(x$gene,rep("None", 100-length(x$gene)))))
colnames(sample_df) <- levels(topDiffCluster$cluster)
sample_df


# In[63]:


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


# In[65]:


write.csv(sample_df,paste0(OS_path_outputs,"AMS_Top100HVGs.csv"),row.names = FALSE)


# ## 8. Cell type annotation

# ### 8.1 mid_KW

# In[32]:


mid <- readRDS("~/Desktop/XMAS_analysis/SL_midbrain_2.rds")


# In[449]:


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


# In[450]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="Cell_type_SL", reduction = "umap_pca_RNA", pt.size=1, label=TRUE) + NoAxes() + ggtitle("RNA")


# In[451]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="orig.ident", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE, cols=orig.ident_colors ) + NoAxes() + ggtitle("Original Identity")
DimPlot(AMS, group.by="Cell_type_SL", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("RNA")


# In[452]:


table(AMS$Cell_type_SL)


# In[453]:


table(AMS$orig.ident, AMS$Cell_type_SL)


# In[454]:


colors <- hue_pal(h.start = 30)(36)

expression <- as.data.frame(as.matrix(AMS@meta.data[,colnames(AMS@meta.data)[colnames(AMS@meta.data) %like% "Cell_type_SL_prediction.score."]]))
colnames(expression) <- gsub("Cell_type_SL_prediction.score.", "", colnames(expression))

expression <- expression[,c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")]
    expression <-cbind(expression,Clusters=AMS$seurat_clusters_pca_RNA)
expression_melt <- reshape2::melt(expression,id=c("Clusters"))
colnames(expression_melt)[3] <- "prediction.score"

expression_melt <- expression_melt %>% mutate(Clusters=factor(Clusters,levels=levels(AMS$seurat_clusters_pca_RNA)))


# In[455]:


options(repr.plot.width=15, repr.plot.height=15)
p <- ggplot(expression_melt, aes(x=Clusters, y=prediction.score,fill=Clusters))
p +geom_violin(scale="width") + geom_boxplot(width=0.1,outlier.shape = NA,position=position_dodge(1),fill="white")+
  theme_bw()+scale_fill_manual(values=colors)+ylim(0,1)+  
theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="right") +facet_grid(variable ~ .,scales="free") +theme(strip.text.y = element_text(angle = 0))


# In[98]:


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

percentages <- tapply(AMS$Cell_type_SL, AMS@meta.data[["seurat_clusters_pca_RNA"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

A_labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[99]:


A_labels #percentage > 3%


# In[462]:


options(repr.plot.width=17, repr.plot.height=5)

n_colors <- 15
color_palette <- viridis::viridis(n_colors, option = "D")

desired_order <- c("794999", "795008", "795011")

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




# In[ ]:





# In[463]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL' == 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type_SL' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# In[468]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$orig.ident, AMS@meta.data$Cell_type_SL))
mat <- mat[,c("DA", "DA0", "Gaba", "GabaNb", "NbM",  "NbML1", "NProg", "OMTN", "OPC", "ProgFP", "Rgl1", "Rgl2", "Rgl3","RN","Sert")]
circos.clear()
par(cex = 1)
chordDiagram(mat, big.gap = 20, small.gap = 2, order = c(colnames(mat), rownames(mat)), grid.col = orig.ident_colors)  


# In[487]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$seurat_clusters_pca_RNA, AMS@meta.data$Cell_type_SL))
mat <- mat[,c("DA", "DA0", "Gaba", "GabaNb", "NbM",  "NbML1", "NProg", "OMTN", "OPC", "ProgFP", "Rgl1", "Rgl2", "Rgl3","RN","Sert")]
circos.clear()
par(cex = 1)
chordDiagram(mat, big.gap = 20, small.gap = 2, order = c(colnames(mat), rownames(mat)))  


# In[469]:


options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(AMS, slot = 'scale.data', group.by = "Cell_type_SL", features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"),ncol=2)


# In[470]:


saveRDS(AMS,paste0(OS_path_outputs, "04_AMS_AFTER_LT.rds"))


# ### 8.2 Developing brain_SL

# In[112]:


mid <- readRDS("~/Desktop/R_data/SL_midbrain_1.rds")


# In[19]:


table(mid$cell.class)


# In[20]:


mid <- SCTransform(mid)
DefaultAssay(AMS) <- "RNA"

Am_anchor  <- FindTransferAnchors(reference = mid, query = AMS, 
                                  normalization.method = 'SCT', recompute.residuals = TRUE,
                                  reduction = 'cca', dims = 1:10)
Am_predictions <- TransferData(anchorset = Am_anchor, refdata = mid$cell.class, dims = 1:10,
                              weight.reduction = 'cca')
AMS <- AddMetaData(AMS, Am_predictions$predicted.id, col.name = 'Cell_type')
for (prediction_score in colnames(Am_predictions)[!colnames(Am_predictions) %in% c("predicted.id", "prediction.score.max")]){
  AMS <- AddMetaData(AMS, Am_predictions[prediction_score], col.name = paste0("Cell_type_",prediction_score))
}


# In[21]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="Cell_type", reduction = "umap_pca_RNA", pt.size=1, label=TRUE) + NoAxes() + ggtitle("RNA")


# In[22]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="orig.ident", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE, cols=orig.ident_colors ) + NoAxes() + ggtitle("Original Identity")
DimPlot(AMS, group.by="Cell_type", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("RNA")


# In[23]:


table(AMS$Cell_type)


# In[24]:


table(AMS$orig.ident, AMS$Cell_type)


# In[32]:


colors <- hue_pal(h.start = 30)(36)

expression <- as.data.frame(as.matrix(AMS@meta.data[,colnames(AMS@meta.data)[colnames(AMS@meta.data) %like% "Cell_type_prediction.score."]]))
colnames(expression) <- gsub("Cell_type_prediction.score.", "", colnames(expression))

expression <- expression[,c("Fibroblast","Glioblast","Neuroblast","Neuron","Oligo","Radial.glia")]
    expression <-cbind(expression,Clusters=AMS$seurat_clusters_pca_RNA)
expression_melt <- reshape2::melt(expression,id=c("Clusters"))
colnames(expression_melt)[3] <- "prediction.score"

expression_melt <- expression_melt %>% mutate(Clusters=factor(Clusters,levels=levels(AMS$seurat_clusters_pca_RNA)))


# In[33]:


options(repr.plot.width=15, repr.plot.height=15)
p <- ggplot(expression_melt, aes(x=Clusters, y=prediction.score,fill=Clusters))
p +geom_violin(scale="width") + geom_boxplot(width=0.1,outlier.shape = NA,position=position_dodge(1),fill="white")+
  theme_bw()+scale_fill_manual(values=colors)+ylim(0,1)+  
theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="right") +facet_grid(variable ~ .,scales="free") +theme(strip.text.y = element_text(angle = 0))


# In[38]:


cell_class <- c("Fibroblast","Glioblast","Neuroblast","Neuron","Oligo","Radial glia")

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

percentages <- tapply(AMS$Cell_type, AMS@meta.data[["orig.ident"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

A_labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[39]:


A_labels #percentage > 3%


# In[102]:


cell_class <- c("Fibroblast","Glioblast","Neuroblast","Neuron","Oligo","Radial glia")

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

percentages <- tapply(AMS$Cell_type, AMS@meta.data[["seurat_clusters_pca_RNA"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

A_labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[103]:


A_labels #percentage > 3%


# In[99]:


A_labels #percentage > 3%


# In[43]:


options(repr.plot.width=17, repr.plot.height=5)

n_colors <- 6
color_palette <- viridis::viridis(n_colors, option = "D")

desired_order <- c("794999", "795008", "795011")

t <- NULL
ggplot(t, aes(x = factor(AMS@meta.data[["orig.ident"]], levels = desired_order), fill = AMS$Cell_type)) +
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




# In[ ]:





# In[49]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type' == 'Neuron',]))) + scale_color_manual(labels = c("Other cell types", "Neuron"), values = c("grey", "darkblue"))
p2 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type' == 'Oligo',]))) + scale_color_manual(labels = c("Other cell types", "Oligo"), values = c("grey", "darkgreen"))
p3 <- DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(AMS@meta.data[AMS@meta.data$'Cell_type' == 'Neuroblast',]))) + scale_color_manual(labels = c("Other cell types", "Neuroblast"), values = c("grey", "red"))
p1 + p2 + p3


# In[45]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$orig.ident, AMS@meta.data$Cell_type))
mat <- mat[,c("Fibroblast","Glioblast","Neuroblast","Neuron","Oligo","Radial glia")]
circos.clear()
par(cex = 1)
chordDiagram(mat, big.gap = 20, small.gap = 2, order = c(colnames(mat), rownames(mat)), grid.col = orig.ident_colors)  


# In[46]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(AMS@meta.data$seurat_clusters_pca_RNA, AMS@meta.data$Cell_type))
mat <- mat[,c("Fibroblast","Glioblast","Neuroblast","Neuron","Oligo","Radial glia")]
circos.clear()
par(cex = 1)
chordDiagram(mat, big.gap = 20, small.gap = 2, order = c(colnames(mat), rownames(mat)))  


# In[47]:


options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(AMS, slot = 'scale.data', group.by = "Cell_type", features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"),ncol=2)


# In[ ]:





# ### 8.3 LaManno

# In[89]:


mid <- readRDS("~/Desktop/XMAS_analysis/LaMannoEmbryo.rds")


# In[90]:


mid <- as.Seurat(mid, data=NULL)


# In[91]:


table(mid$Cell_type)


# In[92]:


mid@assays$RNA <- mid@assays$originalexp 


# In[93]:


mid <- SCTransform(mid)
DefaultAssay(AMS) <- "RNA"

Am_anchor  <- FindTransferAnchors(reference = mid, query = AMS, 
                                  normalization.method = 'SCT', recompute.residuals = TRUE,
                                  reduction = 'cca', dims = 1:10)
Am_predictions <- TransferData(anchorset = Am_anchor, refdata = mid$Cell_type, dims = 1:10,
                              weight.reduction = 'cca')
AMS <- AddMetaData(AMS, Am_predictions$predicted.id, col.name = 'Cell_type_LM')
for (prediction_score in colnames(Am_predictions)[!colnames(Am_predictions) %in% c("predicted.id", "prediction.score.max")]){
  AMS <- AddMetaData(AMS, Am_predictions[prediction_score], col.name = paste0("Cell_type_LM_",prediction_score))
}


# In[94]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="Cell_type_LM", reduction = "umap_pca_RNA", pt.size=1, label=TRUE) + NoAxes() + ggtitle("RNA")


# In[95]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(AMS, group.by="orig.ident", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE, cols=orig.ident_colors ) + NoAxes() + ggtitle("Original Identity")
DimPlot(AMS, group.by="Cell_type_LM", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("RNA")


# In[96]:


table(AMS$Cell_type_LM)


# In[97]:


table(AMS$orig.ident, AMS$Cell_type_LM)


# In[122]:


library(dplyr)

AMS@meta.data <- AMS@meta.data %>%
  select(-starts_with("Cell_type_SL_prediction."))


# ### 8.4 Annotation by clusters

# In[104]:


ml_oc <- c(
    "0_Oligo", 
    "1_Radial glia",
    "2_Neuron",
    "3_Glioblast",
    "4_Glioblast",
    "5_Neuron",
    "6_Neuron",
    "7_Neuron",
    "8_Neuron",
    "9_Neuron",
    "10_Neuron",
    "11_Neuron",
    "12_Neuron",
    "13_Radial glia",
    "14_Oligo",
    "15_Oligo",
    "16_Glioblast", 
    "17_Neuron"
)


# In[105]:


ml_oc <- gsub(".*_", "", ml_oc)


# In[106]:


unique(ml_oc)


# In[107]:


cellType_merged_colors <-  viridis(4)
names(cellType_merged_colors) <- c("Oligo",  "Radial glia", "Neuron",
                                   "Glioblast")


# In[108]:


AMS$cell.class <- ml_oc[AMS@meta.data$seurat_clusters_pca_RNA]
AMS$cell.class <- factor(x = AMS$cell.class, levels = names(cellType_merged_colors))


# In[154]:


table(AMS$cell.class)


# In[111]:


DimPlot(AMS, group.by="cell.class", reduction = "umap_pca_RNA", pt.size=0.5, label=TRUE) + NoAxes() 


# In[136]:


mid <- readRDS("~/Desktop/R_data/SL_midbrain_1.rds")


# In[137]:


library(MetaNeighbor)
library(SummarizedExperiment)


# In[138]:


AMS@meta.data$dataset <- 'graft'
mid$dataset <- "ref"


# In[139]:


Idents(mid)=mid@meta.data$cell.class
Idents(AMS)=AMS@meta.data$cell.class


# In[140]:


mid


# In[141]:


sampled_cells <- mid@meta.data %>%
  rownames() %>%
  sample(size = round(0.6 * ncol(mid)))

# Subset
mid <- subset(mid, cells = sampled_cells)


# In[142]:


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


# In[143]:


Study_ID = rep(c('1', '2'), c(ncol(mid),ncol(AMS)))
Celltype = c(as.character(mid@meta.data$cell.class),as.character(AMS$cell.class))


# In[144]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[145]:


length(unique(mid@meta.data$cell.class))
length(unique(AMS$cell.class))


# In[153]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',4)))
rownames(ann_row)=rownames(celltype_NV)
options(repr.plot.width=12, repr.plot.height=12)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=4,cutree_cols=4,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "average",
         annotation_row = ann_row,
         xtick  = FALSE)


# In[48]:


saveRDS(AMS,paste0(OS_path_outputs, "04_AMS_AFTER_LT.rds"))


# ## 9. Data exploration

# ### 9.1 Subset and reclustering

# In[7]:


neu <- subset(AMS, subset = cell.class=='Neuron')


# In[8]:


neu <- FindVariableFeatures(neu, selection.method = "vst", nfeatures = 2000)


# In[9]:


options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = neu, selection.method = "vst", log = TRUE)


# In[10]:


neu <- ScaleData(neu, features = rownames(neu))


# In[11]:


neu <- RunPCA(neu)


# In[12]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(neu, reduction = "pca", ndims = 50)


# In[13]:


neu <- FindNeighbors(neu, reduction = "pca", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA")


# In[14]:


testing_clusters <- neu
for (clust_res in seq(0.2, 3, by=0.5)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[15]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[16]:


neu <- FindClusters(neu, resolution = 1.0, graph.name="snn_pca_RNA")
neu$seurat_clusters_pca_RNA_2 <- neu$seurat_clusters
neu$seurat_clusters <- NULL


# In[17]:


table(neu$seurat_clusters_pca_RNA_2)


# In[45]:


neu <- RunUMAP(neu, reduction = "pca", dims = 1:40, min.dist = 0.8)


# In[46]:


options(repr.plot.width=7, repr.plot.height=6)
DimPlot(neu, group.by="seurat_clusters_pca_RNA_2", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[47]:


options(repr.plot.width=7, repr.plot.height=6)
DimPlot(neu, group.by="seurat_clusters_pca_RNA", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[ ]:





# ### 9.2 Gene cruating

# In[21]:


options(repr.plot.width=6, repr.plot.height=5)
FeaturePlot(neu, reduction = "umap_pca_RNA", 
            features = c("CALM1"), pt.size=0.05)


# In[48]:


options(repr.plot.width=6, repr.plot.height=5)
FeaturePlot(neu, reduction="umap",
            features = c("CALM1"), pt.size=0.05)


# In[49]:


#VLMC,OPC
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(neu, reduction = "umap", 
            features = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "FBLN1", "PDGFRA", "SOX10", 
                         "OLIG1", "MYRF", "CNP","PLP1"), pt.size=0.05)


# In[50]:


#DA
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(neu, reduction = "umap", 
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1"), pt.size=0.5, alpha=0.7)


# In[239]:


library(scCustomize)


# In[207]:


#DA, slot = "scale.data",
options(repr.plot.width=20, repr.plot.height=25)
VlnPlot(neu, slot = "scale.data", features  =c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", "KCNJ6", "SLC18A2", "SLC6A3","CALB1", "LMO3",'PITX3','ALDH1A1','GAD1','GAD2'),
    pt.size = 0.2, ncol = 2)


# In[590]:


# Top 3 markers per neuron type (modifiable)
neuron_types <- list(
  "Dopaminergic" = c("TH", "SLC6A3", "NR4A2"),
  "Serotonergic" = c("TPH2", "FEV", "SLC6A4"),
  "GABAergic"    = c("GAD1", "GAD2", "SLC32A1"),
  "Glutamatergic" = c("SLC17A6", "LHX9", "SHOX2"),
  "Cholinergic"   = c("CHAT", "SLC18A3", "ACHE"),
  "Noradrenergic" = c("DBH", "SLC6A2", "PNMT")
)


options(repr.plot.width=13, repr.plot.height=4)
# Loop through each neuron type
for (type in names(neuron_types)) {
  markers <- neuron_types[[type]]
  p <- FeaturePlot_scCustom(neu, reduction = "umap", 
            features = markers,  pt.size=0.7,colors_use = viridis_dark_high, num_columns=3)
    
# Combine with main title
    print(wrap_plots(p) + 
            plot_annotation(title = paste(type, "Neurons"),
                          theme = theme(plot.title = element_text(size = 20, face = "bold"))))

    }


# ### 9.3 Marker genes identification

# In[209]:


Idents(neu) <- 'seurat_clusters_pca_RNA_2'


# In[213]:


AMS.markers <- FindAllMarkers(neu, slot = "scale.data", min.pct = 0.05, logfc.threshold = 0.5, only.pos = TRUE)


# In[214]:


AMS.markers_best <- AMS.markers[AMS.markers$p_val_adj < 0.05,]


# In[215]:


AMS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7


# In[217]:


top7


# In[233]:


topDiffCluster <- AMS.markers_best %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


# In[234]:


sample_df <- as.data.frame(lapply(split(topDiffCluster, topDiffCluster$cluster), function(x) c(x$gene,rep("None", 50-length(x$gene)))))
colnames(sample_df) <- levels(topDiffCluster$cluster)
sample_df


# In[235]:


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


# In[236]:


write.csv(sample_df,paste0(OS_path_outputs,"Neu_Top50HVGs.csv"),row.names = FALSE)


# In[ ]:




devtools::install_github(repo = "samuel-marsh/scCustomize", ref = "develop")
# In[25]:


library(scCustomize)


# In[246]:


#DA
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(neu, reduction = "umap", 
            features = c("TH", "NR4A2", "FOXA2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1"), pt.size=0.5, alpha=0.7)


# In[51]:


#DA
options(repr.plot.width=20, repr.plot.height=15)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1"),  pt.size=0.7,colors_use = viridis_dark_high)


# In[59]:


#DA
options(repr.plot.width=20, repr.plot.height=15)
FeaturePlot_scCustom(neu, reduction = "umap", slot = 'scale.data',
            features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),  pt.size=0.7,colors_use = viridis_dark_high) 


# In[29]:


#DA
options(repr.plot.width=20, repr.plot.height=15)
FeaturePlot_scCustom(neu, reduction = "umap_pca_RNA", 
            features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1"),  pt.size=0.5,colors_use = viridis_dark_high)


# In[60]:


#DA
options(repr.plot.width=20, repr.plot.height=15)
FeaturePlot_scCustom(mid, reduction = "UMAP", 
            features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),  pt.size=0.1,colors_use = viridis_dark_high)


# In[55]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("RBFOX3"), pt.size=0.5,colors_use = viridis_dark_high)


# In[418]:


options(repr.plot.width=6, repr.plot.height=5)
FeaturePlot_scCustom(AMS, reduction = "umap_pca_RNA", 
            features = c("RBFOX3"), pt.size=0.2,colors_use = viridis_dark_high)


# In[56]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("DDC"), pt.size=0.3,colors_use = viridis_dark_high)


# In[57]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("MEIS2"), pt.size=0.3,colors_use = viridis_dark_high)


# In[54]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", slot = 'scale.data',
            features = c("TH"), pt.size=0.5,colors_use = viridis_dark_high) + NoAxes()


# In[62]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", slot = 'scale.data',
            features = c("STMN2"), pt.size=0.5,colors_use = viridis_dark_high) + NoAxes()


# In[67]:


#DA, slot = "scale.data",
options(repr.plot.width=15, repr.plot.height=20)
VlnPlot(neu, slot = "scale.data", features  =c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),
    pt.size = 0.1, ncol = 2)


# In[68]:


table(neu$seurat_clusters_pca_RNA_2)


# In[72]:


options(repr.plot.width=20, repr.plot.height=15)
FeaturePlot_scCustom(AMS, reduction = "umap_pca_RNA", 
            features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),  pt.size=0.25, colors_use = viridis_dark_high)


# In[88]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(AMS, reduction = "umap_pca_RNA", pt.size = 0.5, cells.highlight= list(rownames(neu@meta.data[neu@meta.data$'seurat_clusters_pca_RNA_2' %in% c(1,7,10,13),]))) + NoAxes() + scale_color_manual(labels = c("Other cell types", "Dopaminergic Neurons"), values = c("grey", "dark green"))


# In[85]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(neu, reduction = "umap", pt.size = 0.75, cells.highlight= list(rownames(neu@meta.data[neu@meta.data$'seurat_clusters_pca_RNA_2' %in% c(1,7,10,13),]))) + NoAxes() + scale_color_manual(labels = c("Other cell types", "Dopaminergic Neurons"), values = c("grey", "dark green"))


# In[113]:


options(repr.plot.width=10, repr.plot.height=10)
plot <- DotPlot(AMS, group.by = 'seurat_clusters_pca_RNA', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "FBLN1", "PDGFRA", "SOX10", 
                         "OLIG1", "MYRF", "CNP","PLP1"),dot.min = 0.05)
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "New condition", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[102]:


options(repr.plot.width=10, repr.plot.height=10)
plot <- DotPlot(AMS, group.by = 'seurat_clusters_pca_RNA', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),dot.min = 0.05)
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "New condition", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[103]:


options(repr.plot.width=10, repr.plot.height=10)
plot <- DotPlot(neu, group.by = 'seurat_clusters_pca_RNA', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),dot.min = 0.05)
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "New condition", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[107]:


options(repr.plot.width=10, repr.plot.height=10)
plot <- DotPlot(neu, group.by = 'seurat_clusters_pca_RNA_2', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),dot.min = 0.05)
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "New condition", y= "") +  
scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[128]:


library(viridis)
options(repr.plot.width=15, repr.plot.height=10)
scRNAtoolVis::jjDotPlot(object = neu,
          gene = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),
          id = 'seurat_clusters_pca_RNA_2',
          dot.col = rev(viridis(3)),
          dot.min = -3,
          dot.max = 10,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# In[178]:


library(viridis)
options(repr.plot.width=30, repr.plot.height=7)
scRNAtoolVis::jjDotPlot(object = AMS,
          gene = c(
  "SOX2", "VIM", "NES", "PLP1", "EDNRB", "SOX9", "TOP2A", "MKI67", "CENPF", "PTTG1", "PCNA",
  "CORIN", "FOXA1", "FOXA2", "OTX2", "LMX1A", "LMX1B", "ASCL1", "NHLH1", "NEUROD4", "NEUROD1", "NEUROG1",
  "TH", "NR4A2", "PITX3", "EN1", "DDC", "SLC18A2", "SLC6A3", "DBH", "ALDH1A1", "SOX6", "LMO3", 
  "CALB1", "CALB2", "NKX2-1", "BARHL1", "BARHL2", "DBX2", "WNT8B", "IRX3", "IRX5", "PITX2", 
  "PAX6", "LHX1", "LHX2", "LHX9", "LEF1", "TCF7L2", "PAX3", "PAX7", "MEIS2", "BARHL2", "TBR1", 
  "SLC17A6", "GAD1", "GAD2", "DLX1", "DLX5", "DLX6", "ARX", "SOX14", "TLE4", "NKX6-1", "NKX2-2", 
  "POU4F1", "LHX1", "LHX5", "TPH2", "SIM1", "ISL1", "PVALB", "SLC6A4", "SLC17A8", "PDGFRA", 
  "COL1A1", "COL1A2", "LUM", "DCN"),
          id = 'seurat_clusters_pca_RNA',
          dot.col = rev(viridis(3)),
          dot.min = -1,
          dot.max = 9,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# In[111]:


options(repr.plot.width=5, repr.plot.height=4)

FeaturePlot_scCustom(AMS, reduction = "umap_pca_RNA", slot = 'scale.data',
            features = c("mapa"), pt.size=0.25,colors_use = viridis_dark_high) + NoAxes()


# In[116]:


#DA, slot = "scale.data",
options(repr.plot.width=15, repr.plot.height=20)
VlnPlot(AMS, slot = "scale.data", group.by = 'seurat_clusters_pca_RNA', features  =c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "FBLN1", "PDGFRA", "SOX10", 
                         "OLIG1", "MYRF", "CNP","PLP1","MAPT"),
    pt.size = 0.1, ncol = 2)


# In[ ]:





# ### 9.4 Remapping with KW

# In[130]:


mid <- readRDS("~/Desktop/XMAS_analysis/SL_midbrain_2.rds")


# In[131]:


mid <- SCTransform(mid)
DefaultAssay(neu) <- "RNA"

Am_anchor  <- FindTransferAnchors(reference = mid, query = neu, 
                                  normalization.method = 'SCT', recompute.residuals = TRUE,
                                  reduction = 'cca', dims = 1:10)
Am_predictions <- TransferData(anchorset = Am_anchor, refdata = mid$LRprediction_labels, dims = 1:10,
                              weight.reduction = 'cca')
neu <- AddMetaData(neu, Am_predictions$predicted.id, col.name = 'Cell_type_KW')
for (prediction_score in colnames(Am_predictions)[!colnames(Am_predictions) %in% c("predicted.id", "prediction.score.max")]){
  neu <- AddMetaData(neu, Am_predictions[prediction_score], col.name = paste0("Cell_type_KW_",prediction_score))
}


# In[132]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(neu, group.by="Cell_type_KW", reduction = "umap", pt.size=1, label=TRUE) + NoAxes() + ggtitle("RNA")


# In[133]:


library(dplyr)

AMS@meta.data <- AMS@meta.data %>%
  select(-starts_with("Cell_type_KW_prediction."))


# In[138]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(neu, group.by="orig.ident", reduction = "umap", pt.size=1.2, label=FALSE, cols=orig.ident_colors ) + NoAxes() + ggtitle("Original Identity")
DimPlot(neu, group.by="Cell_type_KW", reduction = "umap", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("RNA")


# In[140]:


table(neu$Cell_type_KW)


# In[141]:


table(neu$orig.ident, neu$Cell_type_SL)


# In[143]:


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

percentages <- tapply(neu$Cell_type_KW, neu@meta.data[["seurat_clusters_pca_RNA_2"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

A_labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[144]:


A_labels #percentage > 3%


# In[146]:


options(repr.plot.width=17, repr.plot.height=5)

n_colors <- 15
color_palette <- viridis::viridis(n_colors, option = "D")

desired_order <- c("794999", "795008", "795011")

t <- NULL
ggplot(t, aes(x = factor(neu@meta.data[["orig.ident"]], levels = desired_order), fill = neu$Cell_type_KW)) +
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




# ### 9.5 Data visualization

# In[155]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("NKX2-2"), pt.size=0.5,colors_use = viridis_dark_high)


# In[208]:


options(repr.plot.width=5, repr.plot.height=4)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("SLC17A6"), pt.size=0.5,colors_use = viridis_dark_high)


# In[210]:


options(repr.plot.width=15, repr.plot.height=10)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("SLC17A7", "SLC17A6", "SLC17A8", "GLS", "GRIN1", 
                         "GRIN2A", "GRIN2B", "GRIA1", "GRIA2", "TBR1", "SATB2"), pt.size=0.5,colors_use = viridis_dark_high)


# In[211]:


options(repr.plot.width=15, repr.plot.height=12)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c(  "GAD1", "GAD2",         # Glutamate decarboxylase isoforms
  "SLC32A1",              # VGAT – vesicular GABA transporter
  "DLX1", "DLX2", "DLX5", # Transcription factors for interneuron development
  "LHX6",                 # MGE-derived interneurons
  "SST", "PVALB", "VIP",  # Subtypes of cortical interneurons
  "NPY",                  # Often co-expressed with SST
  "ARX",                  # Interneuron identity gene
  "MAF" # MGE-derived cortical GABAergic neurons
                        ), pt.size=0.5,colors_use = viridis_dark_high)


# In[162]:


options(repr.plot.width=10, repr.plot.height=8)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = c("PITX2","NKX6-1","NKX6-2","NKX2-2"), pt.size=0.5,colors_use = viridis_dark_high)


# In[214]:


a9_markers <- c("TH", "SLC6A3", "SLC18A2", "ALDH1A1", "KCNJ6", "SOX6", "LMX1A", "PITX3", "NR4A2")


# In[215]:


options(repr.plot.width=15, repr.plot.height=10)
FeaturePlot_scCustom(neu, reduction = "umap", 
            features = a9_markers, pt.size=0.5,colors_use = viridis_dark_high)


# In[216]:


library(viridis)
options(repr.plot.width=15, repr.plot.height=10)
scRNAtoolVis::jjDotPlot(object = neu,
          gene = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2", "RBFOX3"),
          id = 'seurat_clusters_pca_RNA_2',
          dot.col = rev(viridis(3)),
          dot.min = -3,
          dot.max = 10,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# ### 9.6 Neu annotation

# In[263]:


ml_oc <- c(
    "0_NProg", 
    "1_DA",
    "2_Glut",
    "3_GABA",
    "4_Glut",
    "5_Glut",
    "6_Glut",
    "7_DA",
    "8_NProg",
    "9_Glut",
    "10_DA",
    "11_Glut",
    "12_Glut",
    "13_DA",
    "14_GABA",
    "15_Glut",
    "16_Glut", 
    "17_Glut",
    "18_GABA"
)


# In[264]:


ml_oc <- gsub(".*_", "", ml_oc)


# In[265]:


unique(ml_oc)


# In[270]:


cellType_merged_colors <-  viridis(4)
names(cellType_merged_colors) <- c("NProg",  "DA", "Glut",
                                   "GABA")


# In[271]:


neu$neu.class <- ml_oc[neu@meta.data$seurat_clusters_pca_RNA_2]
neu$neu.class <- factor(x = neu$neu.class, levels = names(cellType_merged_colors))


# In[272]:


table(neu$neu.class)


# In[278]:


library(viridis)
options(repr.plot.width=12, repr.plot.height=4)
scRNAtoolVis::jjDotPlot(object = neu,
          gene = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2"),
          id = 'neu.class',
          dot.col = rev(viridis(3)),
          dot.min = -3,
          dot.max = 13,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# In[680]:


library(viridis)
options(repr.plot.width=12, repr.plot.height=4)
scRNAtoolVis::jjDotPlot(object = neu,
          gene = c("OTX2", "LMX1A","TH", "DDC", "DRD2", "NR4A2", "EN1", "PBX1", "SLC18A2", "SLC6A3",
                         "SOX6","KCNJ6", "PITX3", "ALDH1A1", "CALB1","GAD1","GAD2","SLC17A6","STMN2"),
          id = 'neu.class',
          dot.col = rev(viridis(3)),
          dot.min = -3,
          dot.max = 13,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# In[ ]:


AMS$cell.class


# In[284]:


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


# In[414]:


table_samples_by_cell_type <- AMS@meta.data %>%
  dplyr::group_by(orig.ident, cell.class) %>%
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(cell.class, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))


# In[415]:


table_samples_by_cell_type


# In[302]:


options(repr.plot.width=4.5, repr.plot.height=7)

temp_labels <- AMS@meta.data %>%
  group_by(orig.ident) %>%
  tally()

table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(sample = factor(orig.ident, levels = levels(AMS@meta.data$orig.ident))) %>%
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

  scale_fill_manual(name = 'Cell type', values = c('#F1BB72', '#68A180', '#D6E7A3', '#57C3F3')) +
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


# In[322]:


desired_order <- c("NProg", "GABA", "DA", "Glut")


# In[323]:


table_samples_by_cell_type <- neu@meta.data %>%
  dplyr::group_by(orig.ident, neu.class) %>%
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(neu.class, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(dplyr::select(., all_of(desired_order)))) %>%
  dplyr::select(c("orig.ident", "total_cell_count", all_of(desired_order)))


# In[324]:


table_samples_by_cell_type


# In[327]:


options(repr.plot.width=4.5, repr.plot.height=7)

temp_labels <- AMS@meta.data %>%
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

  scale_fill_manual(name = 'Neuron type', values = c('#E63863', '#E39A35', '#C1E6F3', '#6778AE')) +
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


# In[343]:


#DA, slot = "scale.data",
options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(neu, slot = "scale.data", group.by = 'neu.class', features  =c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2"),
    pt.size = 0.1, ncol = 4)


# In[349]:


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
options(repr.plot.width=3, repr.plot.height=25)
Idents(neu) <- "neu.class"
StackedVlnPlot(neu, c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2"), pt.size=0.5, cols=my36colors)


# In[329]:


saveRDS(AMS, paste0(OS_path_outputs, "05_AMS_afterAnnotation.rds"))
saveRDS(neu, paste0(OS_path_outputs, "05_neu_afterAnnotation.rds"))


# In[ ]:





# In[ ]:





# In[386]:


DA <- subset(neu, subset = neu.class == "DA")


# In[387]:


DA <- FindVariableFeatures(DA, selection.method = "vst", nfeatures = 2000)


# In[388]:


options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = DA, selection.method = "vst", log = TRUE)


# In[389]:


DA <- ScaleData(DA, features = rownames(DA))


# In[390]:


DA <- RunPCA(DA)


# In[391]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(DA, reduction = "pca", ndims = 50)


# In[392]:


DA <- FindNeighbors(DA, reduction = "pca", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA")


# In[393]:


testing_clusters <- DA
for (clust_res in seq(0.2, 3, by=0.5)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[394]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[395]:


DA <- FindClusters(DA, resolution = 1.0, graph.name="snn_pca_RNA")
DA$seurat_clusters_pca_RNA_3 <- DA$seurat_clusters
DA$seurat_clusters <- NULL


# In[396]:


table(DA$seurat_clusters_pca_RNA_3)


# In[397]:


DA


# In[400]:


DA <- RunUMAP(DA, reduction = "pca", dims = 1:40)


# In[401]:


options(repr.plot.width=7, repr.plot.height=6)
DimPlot(DA, group.by="seurat_clusters_pca_RNA_3", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[403]:


options(repr.plot.width=7, repr.plot.height=6)
DimPlot(DA, group.by="seurat_clusters_pca_RNA_3", reduction="umap", seed=seed, pt.size=1, label=TRUE)


# In[339]:


#DA, slot = "scale.data",
options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(DA, slot = "scale.data", group.by = 'seurat_clusters_pca_RNA_2', features  =c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2"),
    pt.size = 0.1, ncol = 4)


# In[407]:


#DA, slot = "scale.data",
options(repr.plot.width=18, repr.plot.height=15)
VlnPlot(DA, slot = "scale.data", group.by = 'seurat_clusters_pca_RNA_3', features  =c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2"),
    pt.size = 0.1, ncol = 4)


# In[679]:


library(viridis)
options(repr.plot.width=12, repr.plot.height=4)
scRNAtoolVis::jjDotPlot(object = DA,
          gene = c("OTX2", "LMX1A","TH", "DDC", "DRD2", "NR4A2", "EN1", "PBX1", "SLC18A2", "SLC6A3",
                         "SOX6","KCNJ6", "PITX3", "ALDH1A1", "CALB1" ),
          id = 'seurat_clusters_pca_RNA_2',
          dot.col = rev(viridis(3)),
          dot.min = -3,
          dot.max = 13,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# In[406]:


library(viridis)
options(repr.plot.width=12, repr.plot.height=10)
scRNAtoolVis::jjDotPlot(object = DA,
          gene = c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2"),
          id = 'seurat_clusters_pca_RNA_3',
          dot.col = rev(viridis(3)),
          dot.min = -3,
          dot.max = 13,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          ytree = F)


# In[ ]:





# In[350]:


library(Nebulosa)


# In[357]:


options(repr.plot.width=13, repr.plot.height=8)
print(plot_density(neu, c("TH", "DDC", "ALDH1A1", "SLC18A2", "SLC6A3" ), reduction = 'umap', pal='viridis', direction=-1, joint = TRUE) + NoAxes() + NoLegend())


# In[ ]:





# In[ ]:





# In[ ]:





# In[358]:


library(Seurat)
library(SeuratDisk)


# In[359]:


SaveH5Seurat(AMS, filename = paste0(OS_path_outputs,"0325_AMS_annotated.h5Seurat"))


# In[360]:


SeuratDisk::Convert(paste0(OS_path_outputs,"0325_AMS_annotated.h5Seurat"), dest = "h5ad")


# In[361]:


SaveH5Seurat(neu, filename = paste0(OS_path_outputs,"0325_neu_annotated.h5Seurat"))


# In[362]:


SeuratDisk::Convert(paste0(OS_path_outputs,"0325_neu_annotated.h5Seurat"), dest = "h5ad")


# In[ ]:





# In[570]:


DA


# In[464]:


table(AMS$seurat_clusters_pca_RNA)


# In[497]:


geneA_pos_cells <- subset(AMS, subset = seurat_clusters_pca_RNA == 13)


# In[498]:


geneA_pos_cells


# In[499]:


length(colnames(geneA_pos_cells))


# In[508]:


geneB_in_geneA <- sum(GetAssayData(AMS, assay = "RNA")["COL1A1", colnames(geneA_pos_cells)] > 0)
percentage <- (geneB_in_geneA / length(colnames(geneA_pos_cells))) * 100
geneB_in_geneA
length(colnames(geneA_pos_cells))
percentage


# In[609]:


# Define genes
geneA <- "COL1A1"
geneB <- "COL3A1"
geneC <- "PDGFRA"
geneD <- "LUM"

# Get cells expressing each gene
cells_A <- WhichCells(AMS, expression = COL1A1 > 0)
cells_B <- WhichCells(AMS, expression = COL3A1 > 0)
cells_C <- WhichCells(AMS, expression = PDGFRA > 0)
cells_D <- WhichCells(AMS, expression = LUM > 0)

# Find intersection (cells expressing all three)
cells_ABCD <- intersect(intersect(intersect(cells_A, cells_B), cells_C),cells_D)

# Calculate percentage
percentage_ABC <- (length(cells_ABC) / ncol(AMS)) * 100
print(paste0(round(percentage_ABC, 2), "% of cells co-express ", geneA, ", ", geneB, ", ", geneC, ",and ", geneD))


# In[614]:


library(ggVennDiagram)
options(repr.plot.width=6, repr.plot.height=5)
# Create a list of gene-expressing cells
gene_sets <- list(
  A = cells_A,
  B = cells_B,
  C = cells_C,
  D = cells_D
)

# Plot Venn diagram
ggVennDiagram(gene_sets, category.names = c(geneA, geneB, geneC, geneD)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none")


# In[615]:


library(ggVennDiagram)
options(repr.plot.width=6, repr.plot.height=5)
# Create a list of gene-expressing cells
gene_sets <- list(
  A = cells_A,
  B = cells_B,
  C = cells_C
)

# Plot Venn diagram
ggVennDiagram(gene_sets, category.names = c(geneA, geneB, geneC)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none")


# In[467]:


options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(geneA_pos_cells, slot = "data", group.by = 'seurat_clusters_pca_RNA', features  =c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "FBLN1", "PDGFRA", "SOX10", 
                         "OLIG1", "MYRF", "CNP","PLP1"),
    pt.size = 0.1, ncol = 4)


# In[ ]:





# In[441]:


gene_list <- c("DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2")


# In[509]:


results <- data.frame(Gene = character(), Percentage = numeric(), stringsAsFactors = FALSE)
geneA_pos_cells <- WhichCells(DA, expression = TH > 0)

for (gene in gene_list) {
  # Check if the gene exists in the Seurat object to avoid errors
  if (gene %in% rownames(GetAssayData(DA, assay = "RNA"))) {
    geneB_in_geneA <- sum(GetAssayData(DA, assay = "RNA")[gene, geneA_pos_cells] > 0)
    percentage <- (geneB_in_geneA / length(geneA_pos_cells)) * 100
    results <- rbind(results, data.frame(Gene = gene, Percentage = percentage))
  } else {
    warning(paste("Gene", gene, "not found in the Seurat object."))
  }
}

# Print results
print(results)


# In[444]:


options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(DA, slot = "data", group.by = 'neu.class', features  =c("TH", "DDC", "NR4A2", "OTX2", "LMX1A", "EN1", "SOX6","PBX1", 
                         "KCNJ6", "SLC18A2", "SLC6A3","PITX3", "ALDH1A1", "CALB1", "DRD2","GAD1","GAD2","SLC17A6","STMN2"),
    pt.size = 0.1, ncol = 4)


# In[536]:


# Define genes
geneA <- "TH"
geneB <- "SOX6"
geneC <- "NR4A2"

# Get cells expressing each gene
cells_A <- WhichCells(DA, expression = TH > 0)
cells_B <- WhichCells(DA, expression = SOX6 > 0)
cells_C <- WhichCells(DA, expression = NR4A2 > 0)

# Find intersection (cells expressing all three)
cells_ABC <- intersect(intersect(cells_A, cells_B), cells_C)

# Calculate percentage
percentage_ABC <- (length(cells_ABC) / ncol(DA)) * 100
print(paste0(round(percentage_ABC, 2), "% of cells co-express ", geneA, ", ", geneB, ", and ", geneC))


# In[552]:


library(ggVennDiagram)
options(repr.plot.width=4, repr.plot.height=4)
# Create a list of gene-expressing cells
gene_sets <- list(
  A = cells_A,
  B = cells_B,
  C = cells_C
)

# Plot Venn diagram
ggVennDiagram(gene_sets, category.names = c(geneA, geneB, geneC)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none")


# In[543]:


# Define genes
geneA <- "TH"
geneB <- "SOX6"
geneC <- "SLC6A3"

# Get cells expressing each gene
cells_A <- WhichCells(DA, expression = TH > 0)
cells_B <- WhichCells(DA, expression = SOX6 > 0)
cells_C <- WhichCells(DA, expression = SLC6A3 > 0)

# Find intersection (cells expressing all three)
cells_ABC <- intersect(intersect(cells_A, cells_B), cells_C)

# Calculate percentage
percentage_ABC <- (length(cells_ABC) / ncol(DA)) * 100
print(paste0(round(percentage_ABC, 2), "% of cells co-express ", geneA, ", ", geneB, ", and ", geneC))


# In[548]:


options(repr.plot.width=4, repr.plot.height=4)
# Create a list of gene-expressing cells
gene_sets <- list(
  A = cells_A,
  B = cells_B,
  C = cells_C
)

# Plot Venn diagram
ggVennDiagram(gene_sets, category.names = c(geneA, geneB, geneC)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none")


# In[555]:


options(repr.plot.width=4, repr.plot.height=4)
# Create a list of gene-expressing cells
gene_sets <- list(
  A = cells_A,
  B = cells_B,
  C = cells_C
)

# Plot Venn diagram
ggVennDiagram(gene_sets, category.names = c(geneA, geneB, geneC)) +
  scale_fill_viridis(direction = -1, alpha = 0.5) +
  theme(legend.position = "none")


# In[568]:


# Define genes
geneA <- "TH"
geneB <- "SOX6"
geneC <- "CALB1"

# Get cells expressing each gene
cells_A <- WhichCells(DA, expression = TH > 0)
cells_B <- WhichCells(DA, expression = SOX6 > 0)
cells_C <- WhichCells(DA, expression = CALB1 > 0)

# Find intersection (cells expressing all three)
cells_ABC <- intersect(intersect(cells_A, cells_B), cells_C)

# Calculate percentage
percentage_ABC <- (length(cells_ABC) / ncol(DA)) * 100
print(paste0(round(percentage_ABC, 2), "% of cells co-express ", geneA, ", ", geneB, ", and ", geneC))


# In[569]:


options(repr.plot.width=4, repr.plot.height=4)
# Create a list of gene-expressing cells
gene_sets <- list(
  A = cells_A,
  B = cells_B,
  C = cells_C
)

# Plot Venn diagram
ggVennDiagram(gene_sets, category.names = c(geneA, geneB, geneC)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none")

