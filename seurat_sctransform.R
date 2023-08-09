#!/usr/bin/env Rscript

# Packages
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(tibble)
library(glue)

library(ggplot2)
library(plotly)
library(RColorBrewer)

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(SeuratData)
library(SeuratWrappers)


# Input Files
current_path <- file.path("results/day0_day50/seurat")
current_obj <- "iPSC_motor_neuron_day0-day50_time-course_filtered_merged_soupX_seurat_obj_v3.RDS" #v3 Seurat object
merged_filt_obj <- readRDS(file.path(current_path, current_obj))


# SCT
VERSION="v3"
options(Seurat.object.assay.version = VERSION)
Sys.getenv("R_LIBS_USER")
Sys.getenv("TMPDIR")

# remove MT genes
merged_filt_obj <- merged_filt_obj[grep("^MT-", rownames(merged_filt_obj), invert = TRUE), ]

#regularized negative binomial regression to normalize UMI count data
#Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
merged_norm_obj <- SCTransform(merged_filt_obj,
                               assay="RNA",
                               new.assay.name = "SCT",
                               vars.to.regress = "percent.mt",
                               method = "glmGamPoi",
                               return.only.var.genes = TRUE,
                               verbose = TRUE)
# try({
#   merged_norm_obj <- merged_norm_obj %>% 
#     RunPCA(npcs = 50, verbose = TRUE) 
# })
# 
# try({
#   merged_norm_obj <- merged_norm_obj %>% 
#     RunUMAP(reduction = "pca",
#             dims = 1:50, 
#             n.components = 3,
#             n.neighbors = 30,# default 30
#             metric = "cosine",
#             min.dist = 0.2, # default 0.3
#             return.model = TRUE,
#             verbose = T) 
#   
#   merged_norm_obj <- merged_norm_obj %>% 
#     FindNeighbors(reduction = "pca",
#                   dims = 1:50, 
#                   k.param = 20,
#                   nn.method = "annoy",
#                   annoy.metric = "cosine",
#                   n.trees = 100,
#                   verbose = T)
#   
#   merged_norm_obj <- merged_norm_obj %>% 
#     FindClusters(resolution = 1.0, 
#                  algorithm = 2, 
#                  verbose = T) 
#   
#   merged_norm_obj <- merged_norm_obj %>% 
#     FindClusters(resolution = 1.5, 
#                  algorithm = 2, 
#                  verbose = T)
#   
#   merged_norm_obj <- AddMetaData(merged_norm_obj,
#                                  as.data.frame(merged_norm_obj@reductions$umap@cell.embeddings))
#   
# })

# save R object
outrds <- glue::glue("results/day0_day50/seurat/iPSC_motor_neuron_day0-day50_time-course_filtered_merged_soupX_normalized_seurat_obj_{VERSION}.RDS")
saveRDS(merged_norm_obj, outrds)
message(glue::glue("completed saving {basename(outrds)}"))




