#!/usr/bin/env python
# coding: utf-8

# # 0. Preprocessing
#!/bin/bash -l  

#SBATCH -A naiss2025-22-366
#SBATCH -t 1-00:00:00
#SBATCH -J VELO_10X_all
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --no-requeue
#SBATCH --output=VELO_10X.out
#SBATCH --error=VELO_10X.err 
#SBATCH --mail-type=END
#SBATCH --mail-user=guochang.lyu@ki.se
#SBATCH --partition main

## Record the start time
start=`date +%s`
## Record the host being run on
echo "Hostname: $(eval hostname)"

cd /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/
mkdir -p "velocyto_10x" && cd "$_"

name=("2" "5" "6")

cellranger_count_sample_

for index in "${!name[@]}";
do

mkdir -p "velocyto_10x_${name[$index]}" && cd "$_"

  call="velocyto run10x -vv -m /cfs/klemming/projects/supr/gclyu07d/Application/hg38_rmsk.gtf \
        /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/cellranger_count/cellranger_count_sample_${name[$index]} \
        /cfs/klemming/projects/supr/gclyu07d/Application/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"

## Echo the call
echo $call
## Evaluate the call
eval $call

cd ..

done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
#!/bin/bash -l  

#SBATCH -A naiss2025-22-366
#SBATCH -t 4:00:00
#SBATCH -J FILTER_BAM
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --no-requeue
#SBATCH --output=FILTER_BAM.out
#SBATCH --error=FILTER_BAM.err 
#SBATCH --mail-type=END
#SBATCH --mail-user=guochang.lyu@ki.se
#SBATCH --partition main

## Record the start time
start=`date +%s`
## Record the host being run on
echo "Hostname: $(eval hostname)"


cd /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/
mkdir -p "filtered_bam_files" && cd "$_"

name=("2" "5" "6")

for index in "${!name[@]}";
do

  call="subset-bam --bam=/cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/cellranger_count/cellranger_count_sample_${name[$index]}/sample_${name[$index]}/outs/possorted_genome_bam.bam \b
        --cell-barcodes=/cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/filtered_barcode_${name[$index]}.csv \
        --out-bam=/cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/filtered_bam_files/filtered_genome_bam_${name[$index]}.bam \
        --cores=16"


## Echo the call
echo $call
## Evaluate the call
eval $call

done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime#!/bin/bash -l  

#SBATCH -A naiss2025-22-366
#SBATCH -t 12:00:00
#SBATCH -J VELO_10X_filtered
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --no-requeue
#SBATCH --output=VELO_10X_filtered.out
#SBATCH --error=VELO_10X_filtered.err 
#SBATCH --mail-type=END
#SBATCH --mail-user=guochang.lyu@ki.se
#SBATCH --partition main


## Record the start time
start=`date +%s`
## Record the host being run on
echo "Hostname: $(eval hostname)"

cd /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/
mkdir -p "velocyto_10x_filtered" && cd "$_"

name=("2" "5" "6")

for index in "${!name[@]}";
do

  call="velocyto run -v -m /cfs/klemming/projects/supr/gclyu07d/Application/hg38_rmsk.gtf \
        -b /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/filtered_barcode_${name[$index]}.csv \
        -o /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/velocyto_10x_filtered/ \
        /cfs/klemming/projects/supr/gclyu07s/25_AMS_SEQ/filtered_bam_files/filtered_genome_bam_${name[$index]}.bam \
        /cfs/klemming/projects/supr/gclyu07d/Application/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"


## Echo the call
echo $call
## Evaluate the call
eval $call


done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
# In[38]:


BiocManager::install('Seurat')


# In[4]:


library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


# In[1]:


remove.packages('SeuratObject')
remove.packages('Seurat')


# In[3]:


install.packages('~/Downloads/seurat-object-4.1.4.tar.gz')
install.packages('~/Downloads/seurat-4.4.0.tar.gz')


# In[5]:


#Set up global parameters
OS_path <- "/Users/gclyu07/Desktop/R_data/"
OS_path_datasets <- paste0(OS_path, "dataset/")
OS_path_inputs  <- paste0(OS_path, "inputs/")
OS_path_outputs <- paste0(OS_path, "outputs/")

seed  <- 1121
options(repr.plot.width=16, repr.plot.height=12)
options(future.globals.maxSize = 8000 * 1024^2)

options(repr.matrix.max.rows=100, repr.matrix.max.cols=100)


# In[6]:


library(Seurat)
library(SeuratDisk)


# In[7]:


sessionInfo()


# In[8]:


neu <- LoadH5Seurat("~/Desktop/R_data/outputs/0325_neu_annotated.h5Seurat")


# In[9]:


sample.list <- SplitObject(neu, split.by = "orig.ident")


# In[10]:


sample.list

for (index in 1:length(sample.list)){
    t <- list()
    for (colname in colnames(sample.list[[index]])){
        i <- unlist(strsplit(colname, '-'))
        new_colname <- paste0(i[1], '-', 1)
        t[length(t) + 1] <- new_colname
}
    print(length(t))
    list.strings <- paste(t, collapse = "\n")
    writeLines(list.strings, con =paste0(OS_path_outputs,'filtered_barcode_',index,".csv"))
}
# In[11]:


library(rhdf5)

test <- H5Fopen('~/Desktop/R_data/outputs/filtered_genome_bam_2_GZ78W.loom', flags="H5F_ACC_RDONLY")cellid <- h5read(test,"col_attrs/CellID")
gene <- h5read(test,"row_attrs/Gene")
expression <- h5read(test,"matrix")dim(cellid)
dim(gene)
dim(expression)length(unique(cellid))
length(unique(gene))velocyto_filtered_795011 <- SeuratDisk::Connect('~/Desktop/R_data/outputs/filtered_genome_bam_2_2910T.loom', mode = 'r+')
velocyto_filtered_794999 <- SeuratDisk::Connect('~/Desktop/R_data/outputs/filtered_genome_bam_5_05RP5.loom', mode = 'r+')
velocyto_filtered_795008 <- SeuratDisk::Connect('~/Desktop/R_data/outputs/filtered_genome_bam_6_3WCXS.loom', mode = 'r+')velocyto_filtered_795011[["matrix"]][1,]velocyto_filtered_795011[["col_attrs/CellID"]][1]velocyto_filtered_795011[["row_attrs"]]velocyto_filtered_795011[["row_attrs/Gene"]][1]
velocyto_filtered_795011[["row_attrs/Chromosome"]][1]
velocyto_filtered_795011[["row_attrs/Start"]][1]
velocyto_filtered_795011[["row_attrs/End"]][1]
# In[12]:


velocyto_filtered_795011 <- H5Fopen('~/Desktop/R_data/outputs/filtered_genome_bam_2_2910T.loom', flags="H5F_ACC_RDONLY")
velocyto_filtered_794999 <- H5Fopen('~/Desktop/R_data/outputs/filtered_genome_bam_5_05RP5.loom', flags="H5F_ACC_RDONLY")
velocyto_filtered_795008 <- H5Fopen('~/Desktop/R_data/outputs/filtered_genome_bam_6_3WCXS.loom', flags="H5F_ACC_RDONLY")

cellid <- h5read(velocyto_filtered_795011, "col_attrs/CellID")
gene <- h5read(velocyto_filtered_795011, "row_attrs/Gene")
expression <- h5read(velocyto_filtered_795011, "matrix")colnames(expression) <- gene
rownames(expression) <- cellid
expression <- t(expression)
expression <- expression[!duplicated(rownames(expression)),]dim(expression)
# In[71]:


velocyto_filtered <- list(velocyto_filtered_795011,velocyto_filtered_794999,velocyto_filtered_795008)

velocyto_filtered_seurat <- lapply(velocyto_filtered, FUN = function(x){as.Seurat(x)})
# In[72]:


# Apply to each loom file (or HDF5 path) in velocyto_filtered
velocyto_filtered_seurat <- lapply(velocyto_filtered, function(x) {
  
  # Read expression matrix and metadata
  cellid <- h5read(x, "col_attrs/CellID")
  gene <- h5read(x, "row_attrs/Gene")
  expression <- h5read(x, "matrix")
  
  # Assign row and column names
  colnames(expression) <- gene
  rownames(expression) <- cellid
  expression <- t(expression)
  rn <- rownames(expression)
  expression <- expression[!(rn %in% rn[duplicated(rn)]), ]

  # Create Seurat object
  CreateSeuratObject(expression, assay = "unspliced")
})


# In[73]:


velocyto_filtered_seurat[[1]]$orig.ident <- "795011"
velocyto_filtered_seurat[[2]]$orig.ident <- "794999"
velocyto_filtered_seurat[[3]]$orig.ident <- "795008"


# In[74]:


head(colnames(neu[,neu$orig.ident == "795011"]))
head(colnames(neu[,neu$orig.ident == "794999"]))
head(colnames(neu[,neu$orig.ident == "795008"]))


# In[75]:


head(colnames(velocyto_filtered_seurat[[1]]))
head(colnames(velocyto_filtered_seurat[[2]]))
head(colnames(velocyto_filtered_seurat[[3]]))


# In[76]:


library(stringr)


# In[77]:


for (index in 1:length(velocyto_filtered_seurat)){
    t  <- list()
    for (colname in colnames(velocyto_filtered_seurat[[index]])) {
        i <- unlist(strsplit(colname, ':'))
        new_colname <- paste0(str_sub(i[2], end = -2),'-', as.numeric(index))
        t[length(t) + 1] <- new_colname
    }
    velocyto_filtered_seurat[[index]] <- RenameCells(velocyto_filtered_seurat[[index]], new.names = t )
}


# In[78]:


head(colnames(velocyto_filtered_seurat[[1]]))
head(colnames(velocyto_filtered_seurat[[2]]))
head(colnames(velocyto_filtered_seurat[[3]]))


# In[79]:


head(rownames(velocyto_filtered_seurat[[1]]))


# In[80]:


grch38_genes <- grep("^GRCh38-", rownames(velocyto_filtered_seurat[[1]]), value = TRUE)


# In[81]:


velocyto_filtered_seurat <- lapply(velocyto_filtered_seurat, function(x) {
  # Subset object to keep only desired features
  x <- subset(x, features = grch38_genes)
})


# In[82]:


head(rownames(velocyto_filtered_seurat[[2]]))
head(colnames(velocyto_filtered_seurat[[2]]))


# In[83]:


velocyto_filtered_seurat <- lapply(X = velocyto_filtered_seurat, FUN = function(x) {
    DefaultAssay(x) <- "unspliced"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# In[89]:


velocyto_filtered_seurat <- lapply(velocyto_filtered_seurat, function(x) {
  # Clean up gene names
  grch38_genes_2 <- gsub("^GRCh38-", "", rownames(x))
  rownames(x@assays$unspliced@counts) <- grch38_genes_2
  rownames(x@assays$unspliced@data) <- grch38_genes_2 
  return(x)
})


# In[90]:


head(rownames(velocyto_filtered_seurat[[2]]))
head(colnames(velocyto_filtered_seurat[[2]]))


# In[91]:


velocyto_filtered_seurat


# In[92]:


AMS.velocyto <- merge(velocyto_filtered_seurat[[1]], y = velocyto_filtered_seurat[2:3], merge.data = TRUE)


# In[93]:


AMS.velocyto

AMS.velocyto[["unspliced"]] <- JoinLayers(AMS.velocyto[["unspliced"]])AMS.velocyto
# In[94]:


AMS.velocyto <- AMS.velocyto[,colnames(AMS.velocyto) %in% colnames(neu)]


# In[95]:


neu
AMS.velocyto


# In[96]:


AMS.velocyto <- FindVariableFeatures(AMS.velocyto, selection.method = "vst", nfeatures = 2000)


# In[97]:


options(repr.plot.width=10, repr.plot.height=7)
VariableFeaturePlot(object = AMS.velocyto, selection.method = "vst", log = TRUE)


# In[98]:


AMS.velocyto <- ScaleData(AMS.velocyto, features = rownames(AMS.velocyto))


# In[99]:


neu
AMS.velocyto
AMS.velocyto[["unspliced"]]


# In[100]:


table(colnames(AMS.velocyto) %in% colnames(neu))

options(Seurat.object.assay.version = "v3")AMS.velocyto[["unspliced"]] <- as(object = AMS.velocyto[["unspliced"]], Class = "Assay")
# In[101]:


class(AMS.velocyto[["unspliced"]])


# In[102]:


AMS.velocyto


# In[103]:


SaveH5Seurat(AMS.velocyto, filename = "~/Desktop/R_data/outputs/s5_041725_neu_unspliced.h5Seurat")


# In[104]:


SeuratDisk::Convert("~/Desktop/R_data/outputs/s5_041725_neu_unspliced.h5Seurat", dest = "h5ad")


# In[ ]:





# In[105]:


neu@assays[["spliced"]] <- neu@assays[["RNA"]]
neu@assays[["RNA"]] <- NULL


# In[106]:


neu[["unspliced"]] <- AMS.velocyto[["unspliced"]]


# In[ ]:


SaveH5Seurat(neu, filename = "~/Desktop/R_data/outputs/s5_041725_neu_with_unspliced.h5Seurat")


# In[ ]:


SeuratDisk::Convert("~/Desktop/R_data/outputs/s5_041725_neu_with_unspliced.h5Seurat", dest = "h5ad")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




