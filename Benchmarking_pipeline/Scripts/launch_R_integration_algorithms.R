
#Author: Erno Hänninen

#Created: 14.01.2023

#Title: launch_R_Integration_algorithms.R

#Description:
#- Calls the integration algorithms written in R from integration_algorithms.R file. 

#Procedure
#- Takes arguments from nextflow pipeline
#- Reads the integration algorithm parameters from the configuration file
#- Executes the desired integration algorithm by calling it from integration_algorithms.R script
#- stores the output to file
    
#List of non-standard modules:
#- batchelor, seurat, conos

#Usage:
#- This script is launched from the pipeline (processes.nf)


require(batchelor)
require(Seurat)
require(conos)

# Extract the argumates passed from nextflow 
args <- commandArgs(trailingOnly = TRUE)                                                                               
preprocessed_adata_files <- args[[1]]                                                                          
tool <- args[[2]]
output_path <- args[[3]] 
script_dir <- args[[4]]
benchmarking_settings_file <- args[[5]]

# Import the integration algorithms from other script
source(paste0(script_dir,"integration_algorithms.R"))

# Function that calls the integration algoritms
call_algorithms <- function(tool, sobj, batch, hvg, output_path) {

    # Seurat CCA
    if(tool == "seurat_cca"){
        # Subset the features from data based on hvgs
        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
        }else {
            hvg <- rownames(sobj@assays$RNA)
        }
        # Calls seurat_cca function
        out = seurat_cca(sobj, batch, hvg)
        # Save the output
        out[['originalexp']] <- NULL #Delete the originalexp assay, to be sure the integrated assay is read in adata.X
        saveRDS(out, file=paste0(output_path,"seurat_cca_output.rds"))

    }

    # RPCA
    if(tool == "seurat_rpca"){
        # Subset the features from data based on hvgs
        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
        }else {
            hvg <- rownames(sobj@assays$RNA)
        }
        # Calls seuratRPCA function
        out = seuratRPCA(sobj, batch, hvg)
        out[['originalexp']] <- NULL #Delete the originalexp assay, to be sure the integrated assay is read in adata.X
        saveRDS(out, file=paste0(output_path, "seurat_rpca_output.rds"))
    }

    # FastMNN
    if(tool == "fastmnn"){
        # Subset the features from data based on hvgs
        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
            if(length(hvg) > 0){
                sobj <- subset(sobj, features=hvg)
            }
        } 
        # Execute fastmnn 
        fastmnn_sobj=run_fastMNN(sobj, batch)
        saveRDS(fastmnn_sobj, file=paste0(output_path,"fastmnn_output.rds"))

    }

    # Conos
    if(tool == "conos"){
        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
            if(length(hvg) > 0){
                sobj <- subset(sobj, features=hvg)
            }
        } 
        con=runConos(sobj, batch)
        if (file.exists(paste0(output_path,"conos_output.rds"))) {
            #Delete file if it exists
            file.remove(paste0(output_path,"conos_output.rds"))
        }
        meta <- function(sobj) {return(sobj@meta.data)}
        metalist <- lapply(con$samples, meta)
        metaM <- do.call(rbind,unname(metalist))
        saveConosForScanPy(con, output.path=output_path, hdf5_filename="conos_adata.h5", pseudo.pca=TRUE, pca=TRUE, verbose=TRUE, metadata.df=metaM)
    }
}


# Get the path to the preprocessed data
preprocessed_files <- read.csv(file = preprocessed_adata_files)
for(i in 1:nrow(preprocessed_files)){
    if(preprocessed_files[i, ][1] == tool){
        file_path = preprocessed_files[i, ][[2]]
    }
}

# Extract the batch and hvg information from the benchmarking_settings.csv file
benchmarking_settings <- read.csv(file = benchmarking_settings_file)
for(i in 1:nrow(benchmarking_settings)){
    if(benchmarking_settings[i, ][1] == tool){
        batch = benchmarking_settings[i, ][[2]]
        hvg = benchmarking_settings[i, ][[3]]
    }
}

# Read the hvgs from file
if(hvg != "None"){
    hvg = paste0(file_path, "_hvg.RDS")
}

# Read the seurat object
sobj = readRDS(file_path)

# call function that executes the integration algorithm
call_algorithms(tool, sobj, batch, hvg, output_path)