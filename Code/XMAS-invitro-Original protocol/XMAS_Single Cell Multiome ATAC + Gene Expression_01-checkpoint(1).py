#!/usr/bin/env python
# coding: utf-8

# # XMAS_Single Cell Multiome ATAC + Gene Expression 01 
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

# In[4]:


#Set up global parameters
OS_path <- "/Users/gclyu07/Desktop/XMAS_analysis/"
macs2_path = "/Users/gclyu07/anaconda3/envs/SC_v4/bin/macs2"
amulet_path = "/Users/gclyu07/Desktop/XMAS_analysis/packages/AMULET-v1.1/"

OS_path <- '~/Desktop/XMAS_analysis/'
OS_path_datasets <- paste0(OS_path, "dataset/")
OS_path_inputs  <- paste0(OS_path, "inputs/")
OS_path_outputs <- paste0(OS_path, "outputs_no_d0/")

seed  <- 921121
options(repr.plot.width=16, repr.plot.height=12)
options(future.globals.maxSize = 8000 * 1024^2)

options(repr.matrix.max.rows=100, repr.matrix.max.cols=100)


# In[5]:


# set up ATAC annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose =FALSE)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "hg38"


# In[22]:


orig.ident_colors <- c("#F0AD4E", "#D9534F", "#428BCA", "#9933CC", "#66CCCC")
names(orig.ident_colors)  <- c("D11", "D16", "D28","D42", "D56")
show_col(orig.ident_colors)


# ## 2. Loading datasets

# ### 2.1 Load Cellranger aggregation
# 

# In[7]:


data_XMAS  <- Read10X(paste0(OS_path_datasets, "XMAS_cellranger_ARC_aggr_02/filtered_feature_bc_matrix"))


# In[8]:


XMAS <- CreateSeuratObject(data_XMAS$`Gene Expression`)


# In[9]:


XMAS[["ATAC"]] <- CreateChromatinAssay(counts = data_XMAS$Peaks, sep = c(":", "-"), 
                                      fragments = paste0(OS_path_datasets, "XMAS_cellranger_ARC_aggr_02/atac_fragments.tsv.gz"), 
                                       annotation = annotations)


# ### 2.3 Add origin of cell identity

# In[10]:


head(colnames(XMAS))


# In[11]:


t  <- list()
for (colname in colnames(XMAS)) {
        i <- unlist(strsplit(colname, '-'))
        new_colname <- paste0(i[1], '-', as.numeric(i[2]) + 1)
        t[length(t) + 1] <- new_colname
    }


# In[12]:


XMAS <- RenameCells(XMAS, new.names = t )


# In[13]:


head(colnames(XMAS))


# In[10]:


name_XMAS <- seq(1,5,1)
names(name_XMAS) <- c("D11", "D16", "D28","D42", "D56")


# In[11]:


name_XMAS


# In[12]:


XMAS$orig.ident <- unlist(lapply(strsplit(rownames(XMAS@meta.data),"-"), 
                                 function(x) names(name_XMAS[name_XMAS == x[2]])))


# In[13]:


XMAS$orig.ident <- factor(x = XMAS$orig.ident, levels = names(orig.ident_colors))


# In[14]:


table(XMAS$orig.ident)


# ### 2.4 Check number of fragments in peaks sample by sample

# In[15]:


# Check number of fragments in peak
total_fragments_D11 <- CountFragments(paste0(OS_path_datasets, "XMAS_D11/atac_fragments.tsv.gz"))
total_fragments_D16 <- CountFragments(paste0(OS_path_datasets, "XMAS_D16/atac_fragments.tsv.gz"))
total_fragments_D28 <- CountFragments(paste0(OS_path_datasets, "XMAS_D28/atac_fragments.tsv.gz"))
total_fragments_D42 <- CountFragments(paste0(OS_path_datasets, "XMAS_D42/atac_fragments.tsv.gz"))
total_fragments_D56 <- CountFragments(paste0(OS_path_datasets, "XMAS_D56/atac_fragments.tsv.gz"))


# In[16]:


rownames(total_fragments_D11) <-paste0(substr(total_fragments_D11$CB,1,nchar(total_fragments_D11$CB)-2), "-1")
rownames(total_fragments_D16) <- paste0(substr(total_fragments_D16$CB,1,nchar(total_fragments_D16$CB)-2), "-2")
rownames(total_fragments_D28) <- paste0(substr(total_fragments_D28$CB,1,nchar(total_fragments_D28$CB)-2), "-3")
rownames(total_fragments_D42) <- paste0(substr(total_fragments_D42$CB,1,nchar(total_fragments_D42$CB)-2), "-4")
rownames(total_fragments_D56) <- paste0(substr(total_fragments_D56$CB,1,nchar(total_fragments_D56$CB)-2), "-5")


# In[17]:


total_fragments_XMAS <- rbind(total_fragments_D11,total_fragments_D16,
                              total_fragments_D28,total_fragments_D42,total_fragments_D56)


# In[18]:


XMAS@meta.data[c("frequency_count", "mononucleosomal", "nucleosome_free", "reads_count")] <- total_fragments_XMAS[colnames(XMAS[["ATAC"]]), c("frequency_count", "mononucleosomal", "nucleosome_free", "reads_count")]


# ## 3. Remove bad cells based on QC metrics
# 

# ### 3.1 QC calculation
# 

# In[19]:


#Add percentage of mt in each cell
DefaultAssay(XMAS) <- "RNA"
XMAS$percent.mt <- PercentageFeatureSet(XMAS, pattern = "^MT-")


# In[20]:


XMAS$log10GenesUMI <- log10(XMAS$nFeature_RNA) / log10(XMAS$nCount_RNA)


# In[21]:


XMAS$mitoRatio <- XMAS@meta.data$percent.mt / 100


# In[22]:


#Add nucleosome signal and TSS enrichment metrics to each cell
DefaultAssay(XMAS) <- "ATAC"
XMAS <- NucleosomeSignal(XMAS)
XMAS <- TSSEnrichment(XMAS)


# In[23]:


saveRDS(XMAS, paste0(OS_path_outputs, "01_XMAS_afterTSS_v2.rds"))


# ### 3.2 QC visualization

# In[24]:


options(repr.plot.width=30, repr.plot.height=5)
for (sample in 1:length(names(name_XMAS))){
    tmp <- subset(XMAS, subset = orig.ident == names(name_XMAS)[sample])
    Idents(object = tmp) <- "orig.ident"
    print(VlnPlot(object = tmp, features = 
                  c("nCount_RNA", "nCount_ATAC", "percent.mt", "TSS.enrichment", "nucleosome_signal", "frequency_count"), 
                  group.by = "orig.ident", cols=orig.ident_colors, ncol = 6, pt.size = 0.1))
}


# In[25]:


options(repr.plot.width=10, repr.plot.height=5)
RidgePlot(object = XMAS, features = 'nCount_RNA', group.by="orig.ident", slot="counts", cols=orig.ident_colors) + xlim(c(0, 20000)) & NoLegend()
RidgePlot(object = XMAS, features = 'nFeature_RNA', group.by="orig.ident", slot="counts", cols=orig.ident_colors) + xlim(c(0, 5000)) & NoLegend()


# In[26]:


options(repr.plot.width=15, repr.plot.height=8)
FeatureScatter(XMAS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by ="orig.ident", cols=orig.ident_colors)


# In[27]:


options(repr.plot.width=15, repr.plot.height=8)
XMAS@meta.data %>% 
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    scale_fill_manual(values=orig.ident_colors) +
    ggtitle("Cell counts")


# In[28]:


options(repr.plot.width=15, repr.plot.height=8)
XMAS@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_density(color="black", alpha = 0.5) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 600) +
    scale_fill_manual(values=orig.ident_colors) 


# In[29]:


options(repr.plot.width=20, repr.plot.height=15)
XMAS@meta.data %>% 
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


# In[30]:


options(repr.plot.width=15, repr.plot.height=8)
XMAS@meta.data %>%
    ggplot(aes(x=log10GenesUMI, color = orig.ident, fill=orig.ident)) +
    geom_density(color='black', alpha = 0.5) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    scale_fill_manual(values=orig.ident_colors)


# ### 3.3 QC cutoffs

# In[31]:


QC_nCount_RNA_max <- c(100000,200000,25000,20000,50000)
QC_nCount_RNA_min <- c(600,600,600,600,600)
QC_nCount_ATAC_max <- c(100000,500000,300000,120000,100000)
QC_nCount_ATAC_min <- c(1000,1000,1000,1000,1000)
QC_nFeature_RNA_min <- c(250,250,250,250,250)
QC_nucleosome_signal_max <- c(1.5,1.5,1.5,1.5,1.5)
QC_TSS_enrichment_min <- c(2,2,2,2,2)
QC_percentMT_max <- c(10,10,15,15,25)


# In[32]:


XMAS
table(XMAS$orig.ident)


# In[33]:


# filter out low quality cells
good_cells <- c()
for (sample in 1:length(names(name_XMAS))){
    tmp <- subset(XMAS, subset = orig.ident == names(name_XMAS)[sample])
    Idents(object = tmp) <- "orig.ident"
    tmp <- subset(
      x = tmp,
      subset = nCount_ATAC < QC_nCount_ATAC_max[sample] &
        nCount_RNA < QC_nCount_RNA_max[sample] &
        nCount_ATAC > QC_nCount_ATAC_min[sample] &
        nCount_RNA > QC_nCount_RNA_min[sample] &
        nFeature_RNA > QC_nFeature_RNA_min[sample] &
        nucleosome_signal < QC_nucleosome_signal_max[sample] &
        TSS.enrichment > QC_TSS_enrichment_min[sample] &
        percent.mt < QC_percentMT_max[sample]
    )
    good_cells <- c(good_cells, colnames(tmp))
}


# In[34]:


XMAS <- subset(XMAS, cells = good_cells)


# In[35]:


XMAS
table(XMAS$orig.ident)


# In[36]:


saveRDS(XMAS, paste0(OS_path_outputs, "02_XMAS_afterQC_v2.rds"))


# ## 4.Peak calling

# In[37]:


# call peaks using MACS2
peaks <- CallPeaks(XMAS, macs2.path = macs2_path)


# In[38]:


# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(fragments = Fragments(XMAS), features = peaks, cells = colnames(XMAS))


# In[39]:


# create a new assay using the MACS2 peak set and add it to the Seurat object
XMAS[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts, sep = c(":", "-"), 
                    fragments = paste0(OS_path_datasets, "XMAS_cellranger_ARC_aggr_02/atac_fragments.tsv.gz"), 
                    annotation = annotations)


# In[40]:


XMAS <- FRiP(XMAS, assay = 'peaks', total.fragments = 'frequency_count')


# In[41]:


#Remove ATAC assays, saving space and time, we already have Peaks assays
DefaultAssay(XMAS) <- "RNA"
XMAS[['ATAC']] <- NULL


# In[42]:


saveRDS(XMAS, paste0(OS_path_outputs, "03_XMAS_afterPeakCalling_v2.rds"))


# ## 5. Gene accessibility

# In[43]:


DefaultAssay(XMAS) <- "peaks"


# In[44]:


gene.activities <- GeneActivity(XMAS, assay = 'peaks', 
                                extend.upstream = 500, extend.downstream = 0, biotypes = NULL)


# In[45]:


XMAS[['GenePromAcc']] <- CreateAssayObject(counts = gene.activities)


# In[46]:


DefaultAssay(XMAS) <- "GenePromAcc"


# In[47]:


XMAS <- NormalizeData(XMAS, assay = 'GenePromAcc', normalization.method = 'LogNormalize', 
                            scale.factor = median(XMAS$nCount_RNA))


# In[49]:


saveRDS(XMAS, paste0(OS_path_outputs, "04_XMAS_afterGenePromAcc_v2.rds"))


# ## 6. Doublets detection
# 

# ### 6.1 Doublets in RNA using DoubletFinder
# 

# In[50]:


tmp_doublets.list <- SplitObject(XMAS, split.by = "orig.ident")


# In[51]:


tmp_doublets.list


# In[52]:


tmp_doublets.list <- lapply(X = tmp_doublets.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
    x <- ScaleData(x, features = rownames(x))
    x <- RunPCA(x)
    x <- FindNeighbors(x, reduction = "pca", dims = 1:15)
    x <- FindClusters(x, resolution = 1)
    x <- RunUMAP(x, reduction = "pca", dims = 1:15)
})


# In[53]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[1]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[1]))
print(bcmvn)


# In[54]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[2]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[2]))
print(bcmvn)


# In[55]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[3]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[3]))
print(bcmvn)


# In[56]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[4]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[4]))
print(bcmvn)


# In[57]:


#Find best pk for doublet finding
sweep.res <- paramSweep_v3(tmp_doublets.list[[5]], PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(paste("Sample : ", names(tmp_doublets.list)[5]))
print(bcmvn)


# In[58]:


best_pk.v <- c(0.2, 0.21, 0.25, 0.16, 0.01, 0.02)
names(best_pk.v) <- names(tmp_doublets.list)


# In[59]:


for (sample in 1:length(tmp_doublets.list)) {

    homotypic.prop <- modelHomotypic(tmp_doublets.list[[sample]]@meta.data$seurat_clusters)
    nExp_poi <- round(0.15*nrow(tmp_doublets.list[[sample]]@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    pN_value <- 0.25

    tmp_doublets.list[[sample]] <- doubletFinder_v3(tmp_doublets.list[[sample]], PCs = 1:15, pN = pN_value, pK = best_pk.v[sample], nExp = nExp_poi.adj, reuse.pANN = FALSE, sct=FALSE)
}


# In[60]:


DF_df <- c()
for (sample in 1:length(tmp_doublets.list)) {
    DF_tmp <- tmp_doublets.list[[sample]]@meta.data[,(ncol(tmp_doublets.list[[sample]]@meta.data)-1):ncol(tmp_doublets.list[[sample]]@meta.data)]
    colnames(DF_tmp) <- c("pANN", "DF.classifications")
    DF_df <- rbind(DF_df, DF_tmp)
}


# In[61]:


XMAS@meta.data[c("pANN", "DF.classifications")] <- DF_df[colnames(XMAS[["RNA"]]), c("pANN", "DF.classifications")]


# In[62]:


tmp_doublets.list <- NULL


# In[63]:


table(XMAS$DF.classifications)


# ### 6.2 Doublets in ATAC using AMULET

# Download RepeatMasker database for human (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
# 
# Create the human autosome file
# 
# Install python pre requierements for AMULET
# 
# Install AMULET
# 
# Modify FragmentFileOverlapCounter.py line 231 to read per_barcode_metrics.csv for Multiomics : https://github.com/UcarLab/AMULET/issues/11

# In[ ]:


get_ipython().run_cell_magic('', '', 'echo "Hostname: $(eval hostname)"\n\nfor count in 0 11 16 28 42 56\n\ndo\n\ncall="../../packages/AMULET-v1.1/AMULET.sh \\\n../../dataset/XMAS_D${count}/atac_fragments.tsv.gz \\\n../../dataset/XMAS_D${count}/per_barcode_metrics.csv \\\n../../packages/AMULET-v1.1/human_autosomes.txt \\\n../../packages/AMULET-v1.1/blacklist_repeats_segdups_rmsk_hg38.bed \\\n./XMAS_D${count} \\\n../../packages/AMULET-v1.1"\n\necho $call\n\neval $call\n\ndone\n\nend=`date +%s`\nruntime=$((end-start))\necho $runtime\n')


# In[64]:


D11_multiplets <- read.table(paste0(OS_path_outputs,
                                   "AMULET_analysis/XMAS_D11/MultipletCellIds_01.txt")) %>% t() %>% as.vector()
D16_multiplets <- read.table(paste0(OS_path_outputs,
                                   "AMULET_analysis/XMAS_D16/MultipletCellIds_01.txt")) %>% t() %>% as.vector()
D28_multiplets <- read.table(paste0(OS_path_outputs,
                                   "AMULET_analysis/XMAS_D28/MultipletCellIds_01.txt")) %>% t() %>% as.vector()
D42_multiplets <- read.table(paste0(OS_path_outputs,
                                   "AMULET_analysis/XMAS_D42/MultipletCellIds_01.txt")) %>% t() %>% as.vector()
D56_multiplets <- read.table(paste0(OS_path_outputs,
                                   "AMULET_analysis/XMAS_D56/MultipletCellIds_01.txt")) %>% t() %>% as.vector()


# In[65]:


#Generate multiplets list
multiplets.list <- list("D11"=D11_multiplets, 
                        "D16"=D16_multiplets, "D28"=D28_multiplets,
                        "D42"=D42_multiplets, "D56"=D56_multiplets)


# In[66]:


multiplets.list


# In[67]:


XMAS$Amulet_classifications <- rep("Singlet", nrow(XMAS@meta.data))


# In[68]:


amulet_doublets <- c()
for (sample in 1:length(multiplets.list)){
    amulet_doublets <- c(amulet_doublets, paste0(substr(multiplets.list[[sample]],1,nchar(multiplets.list[[sample]])-2), "-",sample))
}


# In[69]:


for (i in 1:nrow(XMAS@meta.data)){
    if (rownames(XMAS@meta.data[i,]) %in% amulet_doublets){
        XMAS@meta.data[i,"Amulet_classifications"] <- "Doublet"
    }
}


# In[70]:


table(XMAS$Amulet_classifications)


# In[71]:


saveRDS(XMAS, paste0(OS_path_outputs, "05_XMAS_afterDoubletsDetection.rds"))


# *** Below started from 010924.***

# ## 7. Processing RNA assay
# 

# ### 7.1 Normalization, merging & scaling of RNA assay

# In[35]:


table(XMAS$orig.ident)


# In[99]:


XMAS <- readRDS(paste0(OS_path_outputs, "05_XMAS_afterDoubletsDetection.rds"))


# In[100]:


XMAS.list <- SplitObject(XMAS, split.by = "orig.ident")


# In[101]:


XMAS.list <- lapply(X = XMAS.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# In[102]:


XMAS <- merge(XMAS.list[[1]], y = XMAS.list[2:5], merge.data = TRUE)


# In[103]:


XMAS <- FindVariableFeatures(XMAS, selection.method = "vst", nfeatures = 2000)


# In[104]:


options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = XMAS, selection.method = "vst", log = TRUE)


# In[105]:


XMAS <- ScaleData(XMAS, features = rownames(XMAS))


# In[106]:


saveRDS(XMAS, paste0(OS_path_outputs, "06_XMAS_afterScaling.rds"))


# ### 7.2 Reducing dimension of RNA assay

# In[107]:


DefaultAssay(XMAS) <- "RNA"


# In[108]:


XMAS <- RunPCA(XMAS)


# ### 7.3 Finding clusters in RNA assay
# 

# In[109]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(XMAS, reduction = "pca", ndims = 50)


# In[110]:


XMAS <- FindNeighbors(XMAS, reduction = "pca", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA")


# In[117]:


testing_clusters <- XMAS
for (clust_res in seq(0.4, 3, by=0.3)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[119]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[120]:


XMAS <- FindClusters(XMAS, resolution = 1.3, graph.name="snn_pca_RNA")
XMAS$seurat_clusters_pca_RNA <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# In[121]:


table(XMAS$seurat_clusters_pca_RNA)


# ### 7.4 Cells projection in 2D

# In[122]:


XMAS <- RunUMAP(XMAS, reduction = "pca", reduction.name="umap_pca_RNA", nn.name="snn_pca_RNA", dims = 1:40)


# In[123]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_pca_RNA", reduction="umap_pca_RNA", seed=seed, pt.size=1, label=TRUE)


# In[124]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[125]:


options(repr.plot.width=25, repr.plot.height=15)
DimPlot(XMAS,  group.by="orig.ident", split.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, ncol=3, label=FALSE)


# In[126]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="seurat_clusters_pca_RNA", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed,  pt.size=1.2, label=FALSE) + NoAxes() 


# In[127]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.5, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# In[128]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.2, label=FALSE)


# ### 7.5 Sample proportion
# 

# In[37]:


t <- table(Cluster=XMAS$seurat_clusters_pca_RNA, Batch=XMAS@meta.data[["orig.ident"]])
t <- t[,rev(names(orig.ident_colors))]
t_percent <- round(prop.table(t, 2) * 100 ,digits = 2)/ncol(t)


# In[38]:


options(repr.plot.width=20, repr.plot.height=12)
barplot(t(t_percent), xlab = "Cluster", ylab="Percentage of cells", 
        legend = TRUE, ylim = c(0, round_any(as.integer(max(rowSums(t_percent)))+5, 5, f = ceiling)), col = rev(orig.ident_colors), args.legend = list(bty = "n", x = "top", ncol = 3))


# In[39]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(XMAS@meta.data$orig.ident, XMAS@meta.data$seurat_clusters_pca_RNA))
for (row in 1:nrow(mat)) mat[row,] <- as.integer(mat[row,]/(sum(mat[row,])/min(rowSums(mat))))
circos.clear()
par(cex = 0.8)

chordDiagram(mat, annotationTrack = c("grid", "axis"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))/3), big.gap = 20, small.gap = 2, order = c(colnames(mat), names(orig.ident_colors)), grid.col=rev(orig.ident_colors))
circos.track(track.index = 1, panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0.3, 0.5))}, bg.border = NA)


# ## 8. Processing ATAC assay via Latent Semantic Indexing (LSI)

# ### 8.1 Frequency normalization

# In[129]:


DefaultAssay(XMAS) <- "peaks"
XMAS <- FindTopFeatures(XMAS, assay = "peaks", min.cutoff = 'q5')
XMAS <- RunTFIDF(XMAS, assay="peaks")


# ### 8.2 Reducing dimension of ATAC assay

# In[130]:


XMAS <- RunSVD(XMAS, assay="peaks")


# ### 8.3 Finding clusters in ATAC assay

# In[131]:


options(repr.plot.width=10, repr.plot.height=7)
DepthCor(XMAS)


# In[132]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(XMAS, reduction = "lsi", ndims = 50)


# In[133]:


XMAS <- FindNeighbors(XMAS, reduction = 'lsi', dims = 2:25, graph.name="snn_ATAC")


# In[134]:


testing_clusters <- XMAS
for (clust_res in seq(0.1, 2, by=0.2)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_ATAC")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[135]:


options(repr.plot.width=15, repr.plot.height=10)
clustree(testing_clusters, prefix = "snn_ATAC_res.") + scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[136]:


XMAS <- FindClusters(XMAS, resolution = 1.5, verbose = FALSE, algorithm = 3, graph.name="snn_ATAC")
XMAS$seurat_clusters_ATAC <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# ### 8.4 Cells projection in 2D

# In[137]:


XMAS <- RunUMAP(XMAS, reduction = 'lsi', reduction.name="umap_ATAC", nn.name="snn_ATAC", dims = 2:25)


# In[138]:


table(XMAS$seurat_clusters_ATAC)


# In[139]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_ATAC", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, pt.size=1, label=FALSE)


# In[55]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[140]:


options(repr.plot.width=25, repr.plot.height=15)
DimPlot(XMAS,  group.by="orig.ident", split.by="orig.ident", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, ncol=3, label=FALSE)


# In[141]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="seurat_clusters_ATAC", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed,  pt.size=1.2, label=FALSE) + NoAxes() 


# In[142]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.5, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# ### 8.5 Sample proportion
# 

# In[143]:


t <- table(Cluster=XMAS$seurat_clusters_ATAC, Batch=XMAS@meta.data[["orig.ident"]])
t <- t[,rev(names(orig.ident_colors))]
t_percent <- round(prop.table(t, 2) * 100 ,digits = 2)/ncol(t)


# In[144]:


options(repr.plot.width=20, repr.plot.height=15)
barplot(t(t_percent), xlab = "Cluster", ylab="Percentage of cells", legend = TRUE, ylim = c(0, round_any(as.integer(max(rowSums(t_percent)))+5, 5, f = ceiling)), col = rev(orig.ident_colors), args.legend = list(bty = "n", x = "top", ncol = 3))


# In[145]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(XMAS@meta.data$orig.ident, XMAS@meta.data$seurat_clusters_ATAC))
for (row in 1:nrow(mat)) mat[row,] <- as.integer(mat[row,]/(sum(mat[row,])/min(rowSums(mat))))
circos.clear()
par(cex = 0.8)

chordDiagram(mat, annotationTrack = c("grid", "axis"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))/3), big.gap = 20, small.gap = 2, order = c(colnames(mat), names(orig.ident_colors)), grid.col=rev(orig.ident_colors))
circos.track(track.index = 1, panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0.3, 0.5))}, bg.border = NA)


# ## 9. Joint UMAP visualization

# ### 9.1 Find multi modal neighbors & clusters
# 

# In[146]:


# build a joint neighbor graph using both assays
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


# In[147]:


testing_clusters <- XMAS
for (clust_res in seq(0.1, 2, by=0.2)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "wsnn")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[148]:


options(repr.plot.width=15, repr.plot.height=10)
clustree(testing_clusters, prefix = "wsnn_res.") + scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[149]:


XMAS <- FindClusters(XMAS, resolution = 1.3, verbose = FALSE, graph.name="wsnn")
XMAS$seurat_clusters_BiMod <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# In[150]:


table(XMAS$seurat_clusters_BiMod)


# ### 9.2 Cells projection in 2D
# 

# In[151]:


# build a joint UMAP visualization
XMAS <- RunUMAP(XMAS, nn.name = "weighted.nn", assay = "RNA", reduction.name="umap_BiMod", verbose = TRUE)


# In[152]:


options(repr.plot.width=11, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_BiMod", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[153]:


options(repr.plot.width=11, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_pca_RNA", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[154]:


options(repr.plot.width=11, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_ATAC", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[75]:


XMAS[,XMAS@meta.data$seurat_clusters_BiMod == 19]$nCount_ATAC


# In[155]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_BiMod", shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[156]:


options(repr.plot.width=25, repr.plot.height=15)
DimPlot(XMAS,  group.by="orig.ident", split.by="orig.ident", reduction="umap_BiMod", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=0.5, ncol=3, label=FALSE)


# In[158]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="seurat_clusters_BiMod", reduction="umap_BiMod", 
        shuffle=TRUE, seed=seed,  pt.size=1.2, label=FALSE) + NoAxes() 


# In[159]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_BiMod", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# ### 7.5 Sample proportion

# In[160]:


t <- table(Cluster=XMAS$seurat_clusters_BiMod, Batch=XMAS@meta.data[["orig.ident"]])
t <- t[,rev(names(orig.ident_colors))]
t_percent <- round(prop.table(t, 2) * 100 ,digits = 2)/ncol(t)


# In[161]:


options(repr.plot.width=20, repr.plot.height=15)
barplot(t(t_percent), xlab = "Cluster", ylab="Percentage of cells", legend = TRUE, ylim = c(0, round_any(as.integer(max(rowSums(t_percent)))+5, 5, f = ceiling)), col = rev(orig.ident_colors), args.legend = list(bty = "n", x = "top", ncol = 3))


# In[162]:


options(repr.plot.width=15, repr.plot.height=15)
set.seed(115)
mat = as.matrix(table(XMAS@meta.data$orig.ident, XMAS@meta.data$seurat_clusters_BiMod))
for (row in 1:nrow(mat)) mat[row,] <- as.integer(mat[row,]/(sum(mat[row,])/min(rowSums(mat))))
circos.clear()
par(cex = 0.8)

chordDiagram(mat, annotationTrack = c("grid", "axis"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))/3), big.gap = 20, small.gap = 2, order = c(colnames(mat), names(orig.ident_colors)), grid.col=rev(orig.ident_colors))
circos.track(track.index = 1, panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0.3, 0.5))}, bg.border = NA)


# In[6]:


XMAS <- readRDS(paste0(OS_path_outputs, "07_XMAS_afterBiMod.rds"))


# In[36]:


umap_embeddings <- list(RNA=as.data.frame(XMAS@reductions[['umap_pca_RNA']]@cell.embeddings), 
                        BIMOD=as.data.frame(XMAS@reductions[['umap_BiMod']]@cell.embeddings),
                        ATAC=as.data.frame(XMAS@reductions[['umap_ATAC']]@cell.embeddings))


# In[37]:


for(i in 1:length(names(umap_embeddings))){
  colnames(umap_embeddings[[i]]) <- c("UMAP_1", "UMAP_2")
  umap_embeddings[[i]]$UMAP_1 <- umap_embeddings[[i]]$UMAP_1 + (i-1)*40
  umap_embeddings[[i]]$modality <- names(umap_embeddings)[i]
  umap_embeddings[[i]]$cluster <- XMAS$orig.ident
  umap_embeddings[[i]]$cell_barcode <- rownames(umap_embeddings[[i]])
  umap_embeddings[[i]] <- umap_embeddings[[i]][sample(1:nrow(umap_embeddings[[i]])), ]
}


# In[38]:


umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)
common.cells <- table(umap.embeddings.merge$cell_barcode)
common.cells <- names(common.cells[common.cells==3])
umap.embeddings.merge <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
umap.embeddings.merge$colors <- mapvalues(as.character(umap.embeddings.merge$cluster), from=names(orig.ident_colors), to=orig.ident_colors)


# In[39]:


color_vector <- umap.embeddings.merge$colors
names(color_vector) <- umap.embeddings.merge$cell_barcode


# In[40]:


options(repr.plot.width=30, repr.plot.height=9)
ggplot(data=umap.embeddings.merge,aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(size=0.2, color=umap.embeddings.merge$colors) + 
  geom_line(data=umap.embeddings.merge,aes(group=cell_barcode, color=cell_barcode),alpha=0.1,linewidth=0.1) +
  scale_color_manual(values=color_vector) +
  theme_classic() + NoAxes() + NoLegend()


# In[163]:


saveRDS(XMAS, paste0(OS_path_outputs, "07_XMAS_afterBiMod.rds"))


# ## 10. Cell cycle scoring

# In[164]:


DefaultAssay(XMAS) <- "RNA"


# In[165]:


# CELL CYCLE
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(XMAS, reduction = "umap_BiMod", 
            features = c("PCNA", "TOP2A", "MCM6", "MKI67"), pt.size=0.5)


# In[166]:


cc.genes


# In[167]:


XMAS <- CellCycleScoring(object=XMAS, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)


# In[168]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(XMAS, features = c('S.Score','G2M.Score'), reduction = "umap_BiMod", cols=c("lightgrey", "blue"), pt.size=1, min.cutoff = 0.05)


# In[169]:


cc_colors <- c("#2A75CB", "#F5563D", "#F5AB00")
names(cc_colors)  <- c("S", "G2M", "G1")
show_col(cc_colors)


# In[170]:


XMAS$Phase


# In[171]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(XMAS, features = c('S.Score','G2M.Score'), reduction = "umap_BiMod", cols=c("purple", "yellow"), pt.size=1, min.cutoff = 0.05)


# In[172]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="Phase", reduction="umap_BiMod", shuffle=TRUE, seed=seed, cols=cc_colors, pt.size=1, label=FALSE)


# In[173]:


RidgePlot(XMAS, group.by="Phase", features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)


# In[174]:


RidgePlot(XMAS, group.by="orig.ident", features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)


# In[175]:


saveRDS(XMAS, paste0(OS_path_outputs, "08_XMAS_afterCellCycleScore.rds"))


# In[176]:


XMAS

