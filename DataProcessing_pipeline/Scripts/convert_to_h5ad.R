#!/usr/bin/env python

#Author: Erno Hänninen
#Created: 2023-22-01
#Title: convert_to_h5ad.py

#Description:
#    - loads seurat object to R and converts to h5ad file
#Procedure:
#    - load data to R 
#    - Save it in h5Seurat format and convert the h5Seurat file to h5ad
#List of functions:
#    - convert_to_h5ad

#List of non-standard modules:
#    - Seurat, SeuratDisk


library(Seurat)
library(SeuratDisk)

# Get arguments from nextflow
args <- commandArgs(trailingOnly = TRUE)                                                                                
raw_data <- args[[1]]                                                                          
integrated_data <- args[[2]]

# Function that creates the merged seurat object
merge_timepoint_data <- function(data_path){
    #Load data
    seurat_object <- get(load(data_path))
    #Merge timepoint data 
    merged_object = merge(x = seurat_object$CS13, y=c(seurat_object$CS14, seurat_object$CS15, seurat_object$CS22_hypo,
                                            seurat_object$CS22_2_hypo, seurat_object$GW16_hypo, seurat_object$GW18_hypo,
                                            seurat_object$GW19_hypo, seurat_object$GW20_34_hypo, seurat_object$GW22T_hypo1,
                                            seurat_object$GW25_3V_hypo))}

# Function that saves seurat object and converts it to h5ad format
convert_to_h5ad <- function(seurat_object, filename){
    SaveH5Seurat(seurat_object, filename = paste0(getwd(), filename)) 
    Convert(paste0(getwd(), filename), assay="RNA", dest = "h5ad") 
    #Remove the unnecessary h5Seurat files
    file.remove(paste0(getwd(), filename))
}

# Raw timepoint data
merged_object <- merge_timepoint_data(raw_data) 
print("Merged data: ")
print(merged_object)
convert_to_h5ad(merged_object, "/MergedHerb.h5Seurat")

# Integrated dataset
seurat_object <- get(load(integrated_data)) # Load data
convert_to_h5ad(seurat_object, "/IntegratedHerb.h5Seurat") # Convert the loaded data

