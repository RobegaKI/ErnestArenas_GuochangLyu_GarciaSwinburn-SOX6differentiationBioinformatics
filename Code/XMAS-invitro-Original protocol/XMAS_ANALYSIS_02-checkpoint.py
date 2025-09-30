#!/usr/bin/env python
# coding: utf-8

# XMAS_Single Cell Multiome ATAC + Gene Expression_02

# # Loading environments

# ## Loading R packages

# In[1]:


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


# In[2]:


library(GSVA)
library(MetaNeighbor)
library(scRNAtoolVis)

options(repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("SeuratObject", "4.1.4")
remotes::install_version("Seurat", "4.4.0", upgrade = FALSE) 
# In[3]:


sessionInfo()


# ## Functions definition 

# In[4]:


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


# ## Set up R parameters 
# 

# In[5]:


#Set up global parameters
OS_path <- "/Users/gclyu07/Desktop/XMAS_analysis/"
macs2_path = "/Users/gclyu07/anaconda3/envs/SC_v4/bin/macs2"
amulet_path = "/Users/gclyu07/Desktop/XMAS_analysis/packages/AMULET-v1.1/"

OS_path <- '~/Desktop/XMAS_analysis/'
OS_path_datasets <- paste0(OS_path, "dataset/")
OS_path_inputs  <- paste0(OS_path, "inputs/")
OS_path_outputs <- paste0(OS_path, "outputs_no_d0/")

seed  <- 1121
options(repr.plot.width=16, repr.plot.height=12)
options(future.globals.maxSize = 8000 * 1024^2)

options(repr.matrix.max.rows=100, repr.matrix.max.cols=100)


# In[6]:


orig.ident_colors <- c("#F0AD4E", "#D9534F", "#428BCA", "#9933CC", "#66CCCC")
names(orig.ident_colors)  <- c("D11", "D16", "D28","D42", "D56")
show_col(orig.ident_colors)


# In[7]:


table(XMAS$Amulet_classifications)
table(XMAS$DF.classifications)


# In[6]:


XMAS <- subset(XMAS, subset = Amulet_classifications == 'Singlet' & DF.classifications == 'Singlet')


# In[9]:


table(XMAS$Amulet_classifications)
table(XMAS$DF.classifications)


# In[10]:


table(XMAS$orig.ident)


# ## 2. RNA

# In[7]:


XMAS.list <- SplitObject(XMAS, split.by = "orig.ident")


# In[8]:


XMAS.list <- lapply(X = XMAS.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# In[9]:


XMAS <- merge(XMAS.list[[1]], y = XMAS.list[2:5], merge.data = TRUE)
XMAS <- FindVariableFeatures(XMAS, selection.method = "vst", nfeatures = 2000)
options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = XMAS, selection.method = "vst", log = TRUE)


# In[10]:


XMAS <- ScaleData(XMAS, features = rownames(XMAS))


# In[11]:


DefaultAssay(XMAS) <- "RNA"
XMAS <- RunPCA(XMAS)


# In[12]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(XMAS, reduction = "pca", ndims = 50)


# In[13]:


XMAS <- FindNeighbors(XMAS, reduction = "pca", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA")


# In[15]:


testing_clusters <- XMAS
for (clust_res in seq(0.4, 3, by=0.3)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[16]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[17]:


XMAS <- FindClusters(XMAS, resolution = 1.3, graph.name="snn_pca_RNA")
XMAS$seurat_clusters_pca_RNA <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# In[18]:


table(XMAS$seurat_clusters_pca_RNA)


# In[19]:


XMAS <- RunUMAP(XMAS, reduction = "pca", reduction.name="umap_pca_RNA", nn.name="snn_pca_RNA", dims = 1:40)


# In[20]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[21]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_pca_RNA", reduction="umap_pca_RNA", seed=seed, pt.size=1, label=TRUE)


# In[22]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="seurat_clusters_pca_RNA", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed,  pt.size=1.2, label=FALSE) + NoAxes() 


# In[23]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_pca_RNA", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.5, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# In[43]:


colnames(XMAS@meta.data)


# In[ ]:





# ## 5. Cell cycle

# In[24]:


DefaultAssay(XMAS) <- "RNA"


# In[25]:


# CELL CYCLE
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(XMAS, reduction = "umap_pca_RNA", 
            features = c("PCNA", "TOP2A", "MCM6", "MKI67"), pt.size=0.5)


# In[26]:


XMAS <- CellCycleScoring(object=XMAS, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)


# In[27]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(XMAS, features = c('S.Score','G2M.Score'), reduction = "umap_pca_RNA", cols=c("lightgrey", "blue"), pt.size=1, min.cutoff = 0.05)


# In[28]:


cc_colors <- c("#2A75CB", "#F5563D", "#F5AB00")
names(cc_colors)  <- c("S", "G2M", "G1")
show_col(cc_colors)


# In[29]:


options(repr.plot.width=20, repr.plot.height=7)
FeaturePlot(XMAS, features = c('S.Score','G2M.Score'), reduction = "umap_pca_RNA", cols=c("purple", "yellow"), pt.size=1, min.cutoff = 0.05)


# In[30]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="Phase", reduction="umap_pca_RNA", shuffle=TRUE, seed=seed, cols=cc_colors, pt.size=1, label=FALSE)


# In[31]:


XMAS@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=XMAS$Phase))+
    theme_minimal()


# In[33]:


XMAS <- ScaleData(XMAS, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(XMAS))


# In[36]:


XMAS <- RunPCA(XMAS, npcs = 50, verbose = FALSE)
XMAS <- FindNeighbors(XMAS, dims = 1:30)
XMAS <- FindClusters(XMAS, resolution = 1.3)
XMAS$seurat_clusters_pca_RNA <- XMAS$seurat_clusters


# In[37]:


XMAS <- RunTSNE(XMAS, dims = 1:30)


# In[38]:


options(repr.plot.width=9, repr.plot.height=7)
DimPlot(XMAS, group.by="orig.ident", reduction="umap", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)
DimPlot(XMAS,reduction = "umap",pt.size=1,group.by = "Phase")


# ## 3. ATAC

# In[39]:


DefaultAssay(XMAS) <- "peaks"
XMAS <- FindTopFeatures(XMAS, assay = "peaks", min.cutoff = 'q5')
XMAS <- RunTFIDF(XMAS, assay="peaks")


# In[40]:


XMAS <- RunSVD(XMAS, assay="peaks")


# In[41]:


options(repr.plot.width=10, repr.plot.height=7)
DepthCor(XMAS)


# In[42]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(XMAS, reduction = "lsi", ndims = 50)


# In[43]:


XMAS <- FindNeighbors(XMAS, reduction = 'lsi', dims = 2:25, graph.name="snn_ATAC")


# In[44]:


testing_clusters <- XMAS
for (clust_res in seq(0.1, 2, by=0.2)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_ATAC")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[45]:


options(repr.plot.width=15, repr.plot.height=10)
clustree(testing_clusters, prefix = "snn_ATAC_res.") + scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[46]:


XMAS <- FindClusters(XMAS, resolution = 1.3, verbose = FALSE, algorithm = 3, graph.name="snn_ATAC")
XMAS$seurat_clusters_ATAC <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# In[47]:


XMAS <- RunUMAP(XMAS, reduction = 'lsi', reduction.name="umap_ATAC", nn.name="snn_ATAC", dims = 2:25, 
                min.dist = 0.5, n.neighbors = 75)


# In[48]:


table(XMAS$seurat_clusters_ATAC)


# In[49]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_ATAC", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, pt.size=1, label=FALSE)


# In[50]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[51]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="seurat_clusters_ATAC", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed,  pt.size=1.2, label=FALSE) + NoAxes() 


# In[52]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_ATAC", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.5, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# ## 4. BiMod

# In[53]:


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


# In[54]:


testing_clusters <- XMAS
for (clust_res in seq(0.1, 2, by=0.2)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "wsnn")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[55]:


options(repr.plot.width=15, repr.plot.height=10)
clustree(testing_clusters, prefix = "wsnn_res.") + scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[56]:


XMAS <- FindClusters(XMAS, resolution = 1.3, verbose = FALSE, graph.name="wsnn")
XMAS$seurat_clusters_BiMod <- XMAS$seurat_clusters
XMAS$seurat_clusters <- NULL


# In[57]:


table(XMAS$seurat_clusters_BiMod)


# In[58]:


# build a joint UMAP visualization
XMAS <- RunUMAP(XMAS, nn.name = "weighted.nn", assay = "RNA", reduction.name="umap_BiMod", verbose = TRUE, 
                min.dist = 1, n.neighbors = 75)


# In[59]:


options(repr.plot.width=11, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_BiMod", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[60]:


options(repr.plot.width=11, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_pca_RNA", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[61]:


options(repr.plot.width=11, repr.plot.height=9)
DimPlot(XMAS, group.by="seurat_clusters_ATAC", reduction="umap_BiMod", seed=seed, pt.size=1, label=TRUE)


# In[62]:


options(repr.plot.width=13, repr.plot.height=9)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_BiMod", shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1, label=FALSE)


# In[63]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="seurat_clusters_BiMod", reduction="umap_BiMod", 
        shuffle=TRUE, seed=seed,  pt.size=1.2, label=FALSE) + NoAxes() 


# In[64]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction="umap_BiMod", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# In[65]:


umap_embeddings <- list(RNA=as.data.frame(XMAS@reductions[['umap_pca_RNA']]@cell.embeddings), 
                        BIMOD=as.data.frame(XMAS@reductions[['umap_BiMod']]@cell.embeddings),
                        ATAC=as.data.frame(XMAS@reductions[['umap_ATAC']]@cell.embeddings))


# In[66]:


for(i in 1:length(names(umap_embeddings))){
  colnames(umap_embeddings[[i]]) <- c("UMAP_1", "UMAP_2")
  umap_embeddings[[i]]$UMAP_1 <- umap_embeddings[[i]]$UMAP_1 + (i-1)*40
  umap_embeddings[[i]]$modality <- names(umap_embeddings)[i]
  umap_embeddings[[i]]$cluster <- XMAS$orig.ident
  umap_embeddings[[i]]$cell_barcode <- rownames(umap_embeddings[[i]])
  umap_embeddings[[i]] <- umap_embeddings[[i]][sample(1:nrow(umap_embeddings[[i]])), ]
}


# In[67]:


umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)
common.cells <- table(umap.embeddings.merge$cell_barcode)
common.cells <- names(common.cells[common.cells==3])
umap.embeddings.merge <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
umap.embeddings.merge$colors <- mapvalues(as.character(umap.embeddings.merge$cluster), from=names(orig.ident_colors), to=orig.ident_colors)


# In[68]:


color_vector <- umap.embeddings.merge$colors
names(color_vector) <- umap.embeddings.merge$cell_barcode


# In[69]:


options(repr.plot.width=30, repr.plot.height=9)
ggplot(data=umap.embeddings.merge,aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(size=0.5, color=umap.embeddings.merge$colors) + 
  geom_line(data=umap.embeddings.merge,aes(group=cell_barcode, color=cell_barcode),alpha=0.1,linewidth=0.15) +
  scale_color_manual(values=color_vector) +
  theme_classic() + NoAxes() + NoLegend()


# ## 6. Marker genes

# In[73]:


DefaultAssay(XMAS) <- "RNA"


# In[74]:


#Markers requested by NW
options(repr.plot.width=20, repr.plot.height=6)
FeaturePlot(XMAS, reduction = "umap_BiMod", 
            features = c("CHCHD2","CHCHD10"))


# In[75]:


#Markers requested by NW
options(repr.plot.width=20, repr.plot.height=12)
FeaturePlot(XMAS, reduction = "umap", slot = 'scale.data', min.cutoff = 0,
            features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"))


# ## 7. Label transfer SL

# In[122]:


mid <- readRDS(paste0(OS_path,"SL_midbrain_2.rds"))


# In[123]:


mid <- SCTransform(mid)
DefaultAssay(XMAS) <- "RNA"

Xm_anchor  <- FindTransferAnchors(reference = mid, query = XMAS, 
                                  normalization.method = 'SCT', recompute.residuals = TRUE,
                                  reduction = 'cca', dims = 1:10)
Xm_predictions <- TransferData(anchorset = Xm_anchor, refdata = mid$LRprediction_labels, dims = 1:10,
                              weight.reduction = 'cca')
XMAS <- AddMetaData(XMAS, Xm_predictions$predicted.id, col.name = 'Cell_type_SL')
for (prediction_score in colnames(Xm_predictions)[!colnames(Xm_predictions) %in% c("predicted.id", "prediction.score.max")]){
  XMAS <- AddMetaData(XMAS, Xm_predictions[prediction_score], col.name = paste0("Cell_type_SL_",prediction_score))
}


# In[78]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(XMAS, group.by="Cell_type_SL", reduction = "umap_BiMod", pt.size=1, label=TRUE) + NoAxes() + ggtitle("BiMod")
DimPlot(XMAS, group.by="Cell_type_SL", reduction = "umap_pca_RNA", pt.size=1, label=TRUE) + NoAxes() + ggtitle("RNA")
DimPlot(XMAS, group.by="Cell_type_SL", reduction = "umap_ATAC", pt.size=1, label=TRUE) + NoAxes() + ggtitle("ATAC")


# In[79]:


options(repr.plot.width=9, repr.plot.height=6)
DimPlot(XMAS, group.by="orig.ident", reduction = "umap_BiMod", pt.size=1.2, label=FALSE, cols=orig.ident_colors ) + NoAxes() + ggtitle("Original Identity")
DimPlot(XMAS, group.by="Cell_type_SL", reduction = "umap_BiMod", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("BiMod")
DimPlot(XMAS, group.by="Cell_type_SL", reduction = "umap_pca_RNA", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("RNA")
DimPlot(XMAS, group.by="Cell_type_SL", reduction = "umap_ATAC", pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("ATAC")


# In[80]:


table(XMAS$orig.ident, XMAS$Cell_type_SL)


# In[81]:


colors <- hue_pal(h.start = 30)(36)

expression <- as.data.frame(as.matrix(XMAS@meta.data[,colnames(XMAS@meta.data)[colnames(XMAS@meta.data) %like% "Cell_type_SL_prediction.score."]]))
colnames(expression) <- gsub("Cell_type_SL_prediction.score.", "", colnames(expression))

expression <- expression[,c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")]
expression <-cbind(expression,Clusters=XMAS$seurat_clusters_BiMod)
expression_melt <- reshape2::melt(expression,id=c("Clusters"))
colnames(expression_melt)[3] <- "prediction.score"

expression_melt <- expression_melt %>% mutate(Clusters=factor(Clusters,levels=levels(XMAS$seurat_clusters_BiMod)))


# In[82]:


options(repr.plot.width=15, repr.plot.height=15)
p <- ggplot(expression_melt, aes(x=Clusters, y=prediction.score,fill=Clusters))
p +geom_violin(scale="width") + geom_boxplot(width=0.1,outlier.shape = NA,position=position_dodge(1),fill="white")+
  theme_bw()+scale_fill_manual(values=colors)+ylim(0,1)+  
theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="right") +facet_grid(variable ~ .,scales="free") +theme(strip.text.y = element_text(angle = 0))


# In[83]:


desired_order <- c("D11", "D16", "D28", "D42","D56")
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

percentages <- tapply(XMAS$Cell_type_SL, XMAS@meta.data[["orig.ident"]], calculate_percentages)
unnested_data <- bind_rows(percentages) %>% cbind(names(percentages),.) %>% dplyr::rename('cell_stage' = 'names(percentages)') %>%
  pivot_longer(cols = -cell_stage, names_to = "cell.type", values_to = "percentage")

labels <- as.data.frame(sapply(percentages, function(p) ifelse(p > 3, paste0(round(p, 1), "%"), " ")))


# In[84]:


labels # percentage > 3%


# In[85]:


options(repr.plot.width=17, repr.plot.height=5)

n_colors <- 15
color_palette <- viridis::viridis(n_colors, option = "D")

desired_order <- c("D11", "D16", "D28", "D42", "D56")

t <- NULL
ggplot(t, aes(x = factor(XMAS@meta.data[["orig.ident"]], levels = desired_order), fill = XMAS$Cell_type_SL)) +
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




# In[45]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(XMAS, reduction = "umap_BiMod", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(XMAS, reduction = "umap_BiMod", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' %like% 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(XMAS, reduction = "umap_BiMod", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# In[44]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(XMAS, reduction = "umap", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(XMAS, reduction = "umap", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' %like% 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(XMAS, reduction = "umap", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' == 'Gaba',]))) + scale_color_manual(labels = c("Other cell types", "hGaba"), values = c("grey", "red"))
p1 + p2 + p3


# In[46]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(XMAS, reduction = "umap_ATAC", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(XMAS, reduction = "umap_ATAC", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' %like% 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(XMAS, reduction = "umap_ATAC", pt.size = 0.5, cells.highlight= list(rownames(XMAS@meta.data[XMAS@meta.data$'Cell_type_SL' == 'Gaba',]))) + scale_color_manual(labels = c("Other cell types", "hGaba"), values = c("grey", "red"))
p1 + p2 + p3


# In[91]:


options(repr.plot.width=15, repr.plot.height=15)
VlnPlot(XMAS, slot = 'data', group.by = "Cell_type_SL", features = c("FOXA2", "LMX1A", "EN1", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX", "NR4A2", "TH","PITX3"),ncol=2)


# In[90]:


table(XMAS$Cell_type_SL)


# In[92]:


saveRDS(XMAS,'01_XMAS_cc+BiMod.rds')


# In[ ]:





# ### MetaNeighbor

# In[7]:


XMAS <- readRDS('01_XMAS_cc+BiMod.rds')


# In[8]:


mid <- readRDS('SL_midbrain_2.rds')


# In[9]:


mid$dataset <- "Sten"


# In[10]:


Idents(mid)=mid@meta.data$LRprediction_labels
mid_neu <- subset(mid, LRprediction_labels %in% c('DA','DA0','Gaba','GabaNb','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3'))


# In[11]:


XMAS@meta.data$dataset <- 'EA_proj'


# In[12]:


Idents(XMAS)=XMAS$seurat_clusters_pca_RNA
DefaultAssay(XMAS) <- "RNA"
XMAS_neu <- subset(XMAS, Cell_type_SL %in% c('DA','DA0','Gaba','GabaNb','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3'))


# In[ ]:





# In[13]:


ob.list <- list(mid_neu, XMAS_neu)


# In[14]:


DefaultAssay(XMAS_neu) <- "RNA"


# In[15]:


XMAS_neu


# In[16]:


mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)


# In[17]:


mda.features


# In[18]:


mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)


# In[19]:


DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)

table(Idents(object = mda.integrated),mda.integrated@meta.data$dataset)


# In[20]:


var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)


# In[21]:


Study_ID = rep(c('r', 'e'), c(ncol(mid_neu),ncol(XMAS_neu)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS_neu$seurat_clusters_pca_RNA))


# In[22]:


library(SummarizedExperiment)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[23]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[24]:


library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


# In[30]:


celltype_NV


# In[33]:


gplots::heatmap.2(celltype_NV,
margins=c(11,23),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 11), rep('darkred',23)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)


# In[34]:


library(pheatmap)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',23)))
rownames(ann_row)=rownames(celltype_NV)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[ ]:





# In[55]:


Study_ID = rep(c('r', 'e'), c(ncol(mid_neu),ncol(XMAS_neu)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS_neu$seurat_clusters_BiMod))


# In[56]:


library(SummarizedExperiment)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[57]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[58]:


library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


# In[59]:


celltype_NV


# In[61]:


gplots::heatmap.2(celltype_NV,
margins=c(11,18),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 11), rep('darkred',18)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)


# In[62]:


library(pheatmap)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',18)))
rownames(ann_row)=rownames(celltype_NV)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[ ]:





# In[35]:


Study_ID = rep(c('r', 'e'), c(ncol(mid_neu),ncol(XMAS_neu)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS_neu$Cell_type_SL))


# In[36]:


library(SummarizedExperiment)
dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[37]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[38]:


library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


# In[39]:


gplots::heatmap.2(celltype_NV,
margins=c(15,15),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 11), rep('darkred',11)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)


# In[40]:


library(pheatmap)

ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',11)))
rownames(ann_row)=rownames(celltype_NV)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[ ]:





# ## 9. Dot plots

# In[88]:


DefaultAssay(XMAS) <- "RNA"
Idents(object = XMAS) <- "seurat_clusters_BiMod"

XMAS.markers <- FindAllMarkers(XMAS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
XMAS.markers_best <- XMAS.markers[XMAS.markers$p_val_adj < 0.05,]
XMAS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7

top7


# In[89]:


top <- XMAS.markers_best %>% filter(!grepl("^RP", gene)) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# In[90]:


top


# In[91]:


options(repr.plot.width=30, repr.plot.height=10)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'seurat_clusters_BiMod', scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)',direction = -1) 


# In[92]:


options(repr.plot.width=30, repr.plot.height=10)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'seurat_clusters_BiMod', scale=TRUE,  cols = c("lightgrey", "red"), dot.scale = 10, 
                features = unique(top$gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") 


# In[67]:


options(repr.plot.width=7, repr.plot.height=7.5)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'seurat_clusters_BiMod', dot.min = 0.05,scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 10, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX","PITX3"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") 


# In[70]:


options(repr.plot.width=7, repr.plot.height=5)
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'seurat_clusters_BiMod', dot.min = 0.05,scale.by= "size", scale=TRUE,  cols = c("lightgrey", "blue"), dot.scale = 5, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX","PITX3"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)') 


# In[68]:


options(repr.plot.width=7, repr.plot.height=5.5)

plot <- DotPlot(XMAS, assay = "RNA", group.by = 'Cell_type_SL', dot.min = 0.05,scale.by= "size", cols = c("lightgrey", "blue"), dot.scale = 7.5, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX","PITX3"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)') 


# In[69]:


options(repr.plot.width=7, repr.plot.height=4)
XMAS.mda <- subset(XMAS, subset = Cell_type_SL %in% c('DA','DA0','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3') )
plot <- DotPlot(XMAS.mda, assay = "RNA", group.by = 'Cell_type_SL', dot.min = 0.05,scale.by= "size", cols = c("lightgrey", "blue"), dot.scale = 7.5, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "SHH", "ASCL1", "NEUROG2", "DCX","PITX3"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)') 


# In[95]:


DefaultAssay(XMAS) <- "RNA"
Idents(object = XMAS) <- "Cell_type_SL"

XMAS.markers <- FindAllMarkers(XMAS, slot = "data", min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)
XMAS.markers_best <- XMAS.markers[XMAS.markers$p_val_adj < 0.05,]
XMAS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7


# In[96]:


top7


# In[115]:


top <- XMAS.markers_best %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)


# In[116]:


top


# In[118]:


options(repr.plot.width=20, repr.plot.height=5.5)
gene <- top$gene[!grepl("^RP", top$gene)]
plot <- DotPlot(XMAS, assay = "RNA", group.by = 'Cell_type_SL', dot.min = 0.05,scale.by= "size", cols = c("lightgrey", "blue"), 
                features = unique(gene))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)') 


# ## 10. abv 0.3

# In[6]:


Xm_predictions <- readRDS('~/Desktop/XMAS_analysis/Xm_predictions.rds')


# In[15]:


t <- ifelse(Xm_predictions$prediction.score.max > 0.3,
        Xm_predictions$predicted.id,
            'Unknown')
XMAS <- AddMetaData(XMAS, t, col.name = 'Cell_type_SL_abv0.3')


# In[16]:


non_un <- subset(XMAS, subset = Cell_type_SL_abv0.3 != "Unknown")


# In[17]:


XMAS
table(XMAS$orig.ident)
non_un
table(non_un$orig.ident)


# In[19]:


options(repr.plot.width=12, repr.plot.height=5)

n_colors <- 14
color_palette <- viridis::viridis(n_colors, option = "D")
desired_order <- c("D11", "D16", "D28", "D42", "D56")

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


# In[20]:


desired_order <- c("D11", "D16", "D28", "D42","D56")
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


# In[ ]:





# #### 3.4.1 Re-project in RNA

# In[21]:


table(non_un$orig.ident)


# In[22]:


DefaultAssay(non_un) <- "RNA"


# In[23]:


non_un <- RunPCA(non_un)
options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(non_un, reduction = "pca", ndims = 50)


# In[24]:


non_un <- FindNeighbors(non_un, reduction = "pca", dims = 1:40, prune.SNN = 0, graph.name="snn_pca_RNA_0.3")


# In[25]:


testing_clusters <- non_un
for (clust_res in seq(0.4, 3, by=0.3)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_pca_RNA_0.3")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[26]:


options(repr.plot.width=15, repr.plot.height=15)
clustree(testing_clusters, prefix = "snn_pca_RNA_0.3_res.") +
scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[27]:


non_un <- FindClusters(non_un, resolution = 1.3, graph.name="snn_pca_RNA_0.3")
non_un$seurat_clusters_pca_RNA_0.3 <- non_un$seurat_clusters
non_un$seurat_clusters <- NULL


# In[28]:


table(non_un$seurat_clusters_pca_RNA_0.3)
table(non_un$orig.ident,non_un$seurat_clusters_pca_RNA_0.3)
table(non_un$seurat_clusters_pca_RNA_0.3,non_un$Cell_type_SL_abv0.3)
table(non_un$orig.ident,non_un$Cell_type_SL_abv0.3)


# In[29]:


non_un <- RunUMAP(non_un, reduction = "pca", reduction.name="umap_pca_RNA_0.3", nn.name="snn_pca_RNA_0.3", dims = 1:40,
                 min.dist = 1, n.neighbors = 75)


# In[30]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=TRUE)


# In[31]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="Cell_type_SL_abv0.3", reduction="umap_pca_RNA_0.3", seed=seed, pt.size=1, label=FALSE)


# In[32]:


options(repr.plot.width=8, repr.plot.height=5)

DimPlot(non_un, group.by="Cell_type_SL_abv0.3", reduction="umap_pca_RNA_0.3",  seed=seed, pt.size=1, label=FALSE) + NoAxes() + ggtitle("Cell types")


# In[33]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="orig.ident", reduction="umap_pca_RNA_0.3", cols=orig.ident_colors, seed=seed, pt.size=1, label=FALSE)


# In[34]:


colors <- hue_pal(h.start = 30)(36)

expression <- as.data.frame(as.matrix(non_un@meta.data[,colnames(non_un@meta.data)[colnames(non_un@meta.data) %like% "Cell_type_SL_prediction.score."]]))
colnames(expression) <- gsub("Cell_type_SL_prediction.score.", "", colnames(expression))

expression <- expression[,c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")]
expression <-cbind(expression,Clusters=non_un$seurat_clusters_pca_RNA_0.3)
expression_melt <- reshape2::melt(expression,id=c("Clusters"))
colnames(expression_melt)[3] <- "prediction.score"

expression_melt <- expression_melt %>% mutate(Clusters=factor(Clusters,levels=levels(non_un$seurat_clusters_pca_RNA_0.3)))

options(repr.plot.width=15, repr.plot.height=15)
p <- ggplot(expression_melt, aes(x=Clusters, y=prediction.score,fill=Clusters))
p +geom_violin(scale="width") + geom_boxplot(width=0.1,outlier.shape = NA,position=position_dodge(1),fill="white")+
  theme_bw()+scale_fill_manual(values=colors)+ylim(0,1)+  
theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="right") +facet_grid(variable ~ .,scales="free") +theme(strip.text.y = element_text(angle = 0))


# In[35]:


colors <- orig.ident_colors

expression <- as.data.frame(as.matrix(non_un@meta.data[,colnames(non_un@meta.data)[colnames(non_un@meta.data) %like% "Cell_type_SL_prediction.score."]]))
colnames(expression) <- gsub("Cell_type_SL_prediction.score.", "", colnames(expression))

expression <- expression[,c("DA","DA0","Gaba","GabaNb","NbM","NbML1","NProg","OMTN","ProgBP","ProgFP","Rgl1","Rgl2","Rgl3","RN")]
expression <-cbind(expression,Clusters=non_un$orig.ident)
expression_melt <- reshape2::melt(expression,id=c("Clusters"))
colnames(expression_melt)[3] <- "prediction.score"

expression_melt <- expression_melt %>% mutate(Clusters=factor(Clusters,levels=unique(non_un$orig.ident)))

options(repr.plot.width=8, repr.plot.height=15)
p <- ggplot(expression_melt, aes(x=Clusters, y=prediction.score,fill=Clusters))
p +geom_violin(scale="width") + geom_boxplot(width=0.1,outlier.shape = NA,position=position_dodge(1),fill="white")+
  theme_bw()+scale_fill_manual(values=colors)+ylim(0,1)+  
theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="right") +facet_grid(variable ~ .,scales="free") +theme(strip.text.y = element_text(angle = 0))


# In[36]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(non_un, reduction = "umap_pca_RNA_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(non_un, reduction = "umap_pca_RNA_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' == 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(non_un, reduction = "umap_pca_RNA_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# #### 3.4.2 Re-project in ATAC

# In[37]:


DefaultAssay(non_un) <- "peaks"


# In[38]:


non_un <- RunSVD(non_un, assay="peaks")


# In[39]:


options(repr.plot.width=10, repr.plot.height=7)
DepthCor(non_un)


# In[40]:


options(repr.plot.width=10, repr.plot.height=7)
ElbowPlot(non_un, reduction = "lsi", ndims = 50)


# In[41]:


non_un <- FindNeighbors(non_un, reduction = 'lsi', dims = 2:25, graph.name="snn_ATAC_0.3")


# In[42]:


testing_clusters <- non_un
for (clust_res in seq(0.1, 2, by=0.2)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "snn_ATAC_0.3")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[43]:


options(repr.plot.width=15, repr.plot.height=10)
clustree(testing_clusters, prefix = "snn_ATAC_0.3_res.") + scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[44]:


non_un <- FindClusters(non_un, resolution = 1.1, verbose = FALSE, algorithm = 3, graph.name="snn_ATAC_0.3")
non_un$seurat_clusters_ATAC_0.3 <- non_un$seurat_clusters
non_un$seurat_clusters <- NULL


# In[45]:


non_un <- RunUMAP(non_un, reduction = 'lsi', reduction.name="umap_ATAC_0.3", nn.name="snn_ATAC_0.3", dims = 2:25,
                                  min.dist = 0.5, n.neighbors = 75)


# In[46]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="seurat_clusters_ATAC_0.3", reduction="umap_ATAC_0.3", seed=seed, pt.size=1, label=TRUE)


# In[47]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="Cell_type_SL_abv0.3", reduction="umap_ATAC_0.3", seed=seed, pt.size=1, label=FALSE)


# In[48]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="orig.ident", reduction="umap_ATAC_0.3", cols=orig.ident_colors, seed=seed, pt.size=1, label=FALSE)


# In[49]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(non_un, reduction = "umap_ATAC_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(non_un, reduction = "umap_ATAC_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' == 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(non_un, reduction = "umap_ATAC_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# #### 3.4.3 Re-project in BiMod

# In[50]:


# build a joint neighbor graph using both assays
non_un <- FindMultiModalNeighbors(
  object = non_un,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:25),
  weighted.nn.name = "weighted.nn",
  modality.weight.name = c("RNA.weight.name", "peaks.weight.name"),
  snn.graph.name = "wsnn",
  prune.SNN = 0,
  verbose = TRUE
)


# In[51]:


testing_clusters <- non_un
for (clust_res in seq(0.1, 2, by=0.2)){
  testing_clusters <- FindClusters(testing_clusters, resolution = clust_res, graph.name = "wsnn")
  if (clust_res > 1 & length(unique(Idents(testing_clusters))) == 1){
    break
  }
}


# In[52]:


options(repr.plot.width=15, repr.plot.height=10)
clustree(testing_clusters, prefix = "wsnn_res.") + scale_edge_color_continuous(low = "lightgrey", high = "darkblue")


# In[53]:


non_un <- FindClusters(non_un, resolution = 1.1, verbose = FALSE, graph.name="wsnn")
non_un$seurat_clusters_BiMod_0.3 <- non_un$seurat_clusters
non_un$seurat_clusters <- NULL


# In[54]:


# build a joint UMAP visualization
non_un <- RunUMAP(non_un, nn.name = "weighted.nn", assay = "RNA", reduction.name="umap_BiMod_0.3", verbose = TRUE,
                  min.dist = 1, n.neighbors = 75)


# In[55]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="seurat_clusters_BiMod_0.3", reduction="umap_BiMod_0.3", seed=seed, pt.size=1, label=TRUE)


# In[56]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="Cell_type_SL_abv0.3", reduction="umap_BiMod_0.3", seed=seed, pt.size=1, label=FALSE)


# In[57]:


options(repr.plot.width=8, repr.plot.height=5)
DimPlot(non_un, group.by="orig.ident", reduction="umap_BiMod_0.3", cols=orig.ident_colors, seed=seed, pt.size=1, label=FALSE)


# In[58]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(non_un, group.by="orig.ident", reduction="umap_BiMod_0.3", 
        shuffle=TRUE, seed=seed, cols=orig.ident_colors, pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("Orignal identity")


# In[59]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(non_un, group.by="Cell_type_SL_abv0.3", reduction="umap_BiMod_0.3", 
        shuffle=TRUE, seed=seed, pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("Cell types")


# In[60]:


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(non_un, group.by="seurat_clusters_BiMod_0.3", reduction="umap_BiMod_0.3", 
        shuffle=TRUE, seed=seed, pt.size=1.2, label=FALSE) + NoAxes() + ggtitle("Cell types")


# In[61]:


options(repr.plot.width=21, repr.plot.height=4)
p1 <- DimPlot(non_un, reduction = "umap_BiMod_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' %like% 'DA',]))) + scale_color_manual(labels = c("Other cell types", "hDA"), values = c("grey", "darkblue"))
p2 <- DimPlot(non_un, reduction = "umap_BiMod_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' %like% 'NbM',]))) + scale_color_manual(labels = c("Other cell types", "NbM"), values = c("grey", "darkgreen"))
p3 <- DimPlot(non_un, reduction = "umap_BiMod_0.3", pt.size = 0.5, cells.highlight= list(rownames(non_un@meta.data[non_un@meta.data$'Cell_type_SL' == 'ProgFP',]))) + scale_color_manual(labels = c("Other cell types", "ProgFP"), values = c("grey", "red"))
p1 + p2 + p3


# In[ ]:





# In[62]:


options(repr.plot.width=9, repr.plot.height=5)
non_un.mda <- subset(non_un, subset = Cell_type_SL %in% c('DA','DA0','NbM','NbML1',
                                 'NProg','ProgFP','Rgl1', 'Rgl2', 'Rgl3') )
plot <- DotPlot(non_un.mda, assay = "RNA", group.by = 'Cell_type_SL', dot.min = 0.03,scale.by= "size", cols = c("lightgrey", "blue"), dot.scale = 7.5, 
                features = c("FOXA2", "LMX1A", "EN1", "NR4A2", "TH","DDC", "OTX2", "SOX6", "LMO3","PBX1","CABL1","KCNJ6","SHH", "ASCL1", "NEUROG2", "DCX","SLC18A2"))
plot + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y= "") +   scale_color_viridis_c(name = 'log2 (count + 1)') 


# In[ ]:





# In[63]:


saveRDS(non_un, "02_abv031_XMAS.rds")


# In[ ]:





# In[64]:


Idents(non_un)=non_un$seurat_clusters_BiMod_0.3
DefaultAssay(non_un) <- "RNA"
XMAS_neu <- non_un


# In[65]:


ob.list <- list(mid_neu, XMAS_neu)


# In[66]:


mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)


# In[67]:


mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)


# In[68]:


DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)

table(Idents(object = mda.integrated),mda.integrated@meta.data$dataset)


# In[83]:


saveRDS(mda.integrated, "mda.integrated_XMAS_NEU_ABV031.rds")


# In[69]:


var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)


# In[84]:


Study_ID = rep(c('r', 'e'), c(ncol(mid_neu),ncol(XMAS_neu)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS_neu$seurat_clusters_BiMod_0.3))


# In[85]:


dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[87]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[88]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


# In[89]:


length(unique(mid_neu@meta.data$LRprediction_labels))
length(unique(XMAS_neu$seurat_clusters_BiMod_0.3))


# In[90]:


options(repr.plot.width=15, repr.plot.height=15)


gplots::heatmap.2(celltype_NV,
margins=c(11,17),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 11), rep('darkred',17)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)


# In[91]:


options(repr.plot.width=8, repr.plot.height=8)



ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',17)))
rownames(ann_row)=rownames(celltype_NV)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[92]:


Study_ID = rep(c('r', 'e'), c(ncol(mid_neu),ncol(XMAS_neu)))
Celltype = c(as.character(mid_neu@meta.data$LRprediction_labels),as.character(XMAS_neu$Cell_type_SL_abv0.3))


# In[93]:


dat=SummarizedExperiment(assays=list(counts=combined_mat))


# In[94]:


celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


# In[95]:


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


# In[96]:


length(unique(mid_neu@meta.data$LRprediction_labels))
length(unique(XMAS_neu$Cell_type_SL_abv0.3))


# In[97]:


gplots::heatmap.2(celltype_NV,
margins=c(11,14),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 11), rep('darkred',14)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)


# In[98]:


ann_row=data.frame(dataset=c(rep("Reference", 11), rep('Experiment',14)))
rownames(ann_row)=rownames(celltype_NV)

pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Reference = "darkgreen", Experiment = "darkred"),
         xtick  = FALSE)
         )


# In[ ]:





# ## 11. clusterprofiling

# In[132]:


options(repr.plot.width=24, repr.plot.height=5)
p1 <- DimPlot(non_un, group.by="seurat_clusters_pca_RNA_0.3", reduction="umap_BiMod_0.3", seed=seed, pt.size=1, label=TRUE)
p2 <- DimPlot(non_un, group.by="seurat_clusters_BiMod_0.3", reduction="umap_BiMod_0.3", seed=seed, pt.size=1, label=TRUE)
p3 <- DimPlot(non_un, group.by="orig.ident", reduction="umap_BiMod_0.3", seed=seed, pt.size=1, label=FALSE, cols = orig.ident_colors)
p1+p2+p3


# In[133]:


table(non_un$seurat_clusters_BiMod_0.3, non_un$Cell_type_SL_abv0.3)
table(non_un$seurat_clusters_BiMod_0.3, non_un$orig.ident)
table(non_un$seurat_clusters_pca_RNA_0.3, non_un$orig.ident)


# In[105]:


Idents(object = non_un) <- "seurat_clusters_BiMod_0.3"
cluster.markers <- FindMarkers(non_un, ident.1 = 4, ident.2 = 16, min.pct = 0.25)  


# In[106]:


head(cluster.markers, n = 5) 


# In[107]:


up <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]>=0.25)),])
down <-rownames(cluster.markers[intersect(which(cluster.markers [,1]<0.05),which(cluster.markers [,2]<=(-0.25))),])


# In[108]:


gs <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)


# In[109]:


ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)                    
head(ego.bp)


# In[113]:


options(repr.plot.width=7, repr.plot.height=5)
print(dotplot(ego.bp, showCategory=10,title="Day16 vs.Day28 (Cluster 4 vs. 16) up gene GoTerm"))


# In[114]:


kk <- enrichKEGG(gene= gs$ENTREZID, organism = 'hsa',pvalueCutoff = 0.05)


# In[115]:


print(dotplot(kk, showCategory=10,title="KEGG_biological")) 


# In[116]:


gs <- bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)


# In[117]:


ego.bp <- enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont= "BP", pAdjustMethod = "BH", pvalueCutoff= 0.05, qvalueCutoff= 0.2,
                   readable= TRUE)                    
head(ego.bp)


# In[119]:


print(dotplot(ego.bp, showCategory=10,title="Day16 vs.Day28 (Cluster 4 vs. 16) down gene GoTerm"))


# In[120]:


kk <- enrichKEGG(gene= gs$ENTREZID, organism = 'hsa',pvalueCutoff = 0.05)


# In[121]:


print(dotplot(kk, showCategory=10,title="KEGG_biological")) 


# ### 11.1 gsea analysis

# In[122]:


library(msigdbr)
library(fgsea)


# In[123]:


markers <- FindMarkers(non_un, ident.1 = 4, ident.2 = 16, min.pct = 0.25, logfc.threshold = 0) 
markers$genes <- rownames(markers)
markers %>% arrange(desc(avg_log2FC)) %>% select(genes,avg_log2FC) -> cluster.genes
ranks<- deframe(cluster.genes)


# In[124]:


mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")


# In[125]:


fgsea_sets <- mdb_c2 [grep("^KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)


# In[126]:


ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") 


# In[191]:


plotEnrichment(fgsea_sets[["KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"]],ranks) 
+ labs(title="KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION")


# In[203]:


geneset <- read.gmt("c8.all.v2023.2.Hs.entrez.gmt") 


# In[204]:


markers <- FindMarkers(non_un, ident.1 = 6, ident.2 = c(1,11), min.pct = 0.1, logfc.threshold = 0)
gs <-bitr(rownames(markers), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
markers1<-cbind(markers[gs[,1],],gs)
geneList <- markers1$avg_log2FC
names(geneList) <- markers1$ENTREZID
geneList <- sort(geneList,decreasing = T)
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F)
egmt1<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")


# In[208]:


options(repr.plot.width=15, repr.plot.height=5)

dotplot(egmt1,split=".sign")+facet_grid(~.sign)


# In[502]:


egmt@result$ID[c(1,3,4)]


# In[507]:


upsetplot(egmt)


# In[571]:


options(repr.plot.width=7, repr.plot.height=5)

gseaplot2(egmt,c(1,3,4),color="red",pvalue_table = F,base_size=7.5,ES_geom="line")


# In[429]:


for(i in seq_along(egmt@result$ID)){
print(gseaplot2(egmt, geneSetID = i, title = egmt@result$ID[i]))
}


# In[ ]:





# In[ ]:





# In[ ]:





# In[196]:


geneset <- read.gmt("c2.cp.kegg.v7.5.entrez.gmt") 


# In[475]:


markers <- FindMarkers(non_un, ident.1 = 3, ident.2 = 1, min.pct = 0.1, logfc.threshold = 0)
gs <-bitr(rownames(markers), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
markers1<-cbind(markers[gs[,1],],gs)
geneList <- markers1$avg_log2FC
names(geneList) <- markers1$ENTREZID
geneList <- sort(geneList,decreasing = T)
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F)
egmt1<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")


# In[476]:


options(repr.plot.width=10, repr.plot.height=10)
dotplot(egmt1,split=".sign")+facet_grid(~.sign)


# In[479]:


options(repr.plot.width=5.5, repr.plot.height=5)
gseaplot2(egmt,"MANNO_MIDBRAIN_NEUROTYPES_HDA",color="red",pvalue_table = F,title="MANNO_MIDBRAIN_NEUROTYPES_HDA Cluster 1 vs. 3",base_size=10,ES_geom="line")


# In[ ]:





# ### 11. jjVolcano

# In[135]:


colnames(non_un@meta.data)


# In[136]:


Idents(non_un) <- "seurat_clusters_BiMod_0.3"


# In[152]:


deg <- FindAllMarkers(non_un, min.pct = 0.1, logfc.threshold = 0.5, only.pos = FALSE)
head(deg)


# In[139]:


library(colorspace)
high_contrast_colors <- rainbow_hcl(17)


# In[153]:


deg


# In[157]:


options(repr.plot.width=30, repr.plot.height=10)

jjVolcano(diffData = deg,
          tile.col = high_contrast_colors,
          pvalue.cutoff = 0.01) + NoLegend()


# In[ ]:





# In[ ]:





# In[158]:


Idents(non_un) <- "Cell_type_SL_abv0.3"


# In[159]:


deg <- FindAllMarkers(non_un, min.pct = 0.1, logfc.threshold = 0.5, only.pos = FALSE)
head(deg)


# In[160]:


library(colorspace)
high_contrast_colors <- rainbow_hcl(17)


# In[161]:


deg


# In[166]:


options(repr.plot.width=15, repr.plot.height=10)
colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
            "#1F78B4", "#33A02C")
jjVolcano(diffData = deg,
          tile.col = colors,
          topGeneN = 3) + NoLegend()


# In[ ]:




