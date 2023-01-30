require(batchelor)
require(Seurat)
require(conos)

args <- commandArgs(trailingOnly = TRUE)  
print(args)                                                                              
preprocessed_adata_files <- args[[1]]                                                                          
tool <- args[[2]]
output_path <- args[[3]] 
script_dir <- args[[4]]
benchmarking_settings_file <- args[[5]]

source(paste0(script_dir,"integration_algorithms.R"))


    
call_algorithms <- function(tool, sobj, batch, hvg, output_path) {

    if(tool == "seurat_cca"){

        print(hvg)
        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
        }else {
            hvg <- rownames(sobj@assays$RNA)
        }
        print(hvg)
        out = seurat_cca(sobj, batch, hvg)

        out[['originalexp']] <- NULL #Delete the originalexp assay, to be sure the integrated assay is read in adata.X
        print("seurat cca output:")
        print(out)
        saveRDS(out, file=paste0(output_path,"seurat_cca_output.rds"))

    }

    if(tool == "seurat_rpca"){

        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
        }else {
            hvg <- rownames(sobj@assays$RNA)
        }

        out = seuratRPCA(sobj, batch, hvg)
        out[['originalexp']] <- NULL #Delete the originalexp assay, to be sure the integrated assay is read in adata.X
        print("seurat rpca output:")
        print(out)
        saveRDS(out, file=paste0(output_path, "seurat_rpca_output.rds"))
    }

    if(tool == "fastmnn"){
        if(file.exists(hvg)) {
            print(hvg)
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
            print(length(hvg))
            if(length(hvg) > 0){
                sobj <- subset(sobj, features=hvg)
            }
        } 
        fastmnn_sobj=run_fastMNN(sobj, batch)
        saveRDS(fastmnn_sobj, file=paste0(output_path,"fastmnn_output.rds"))

    }

    if(tool == "conos"){
        if(file.exists(hvg)) {
            hvg<-unlist(readRDS(hvg), use.names=FALSE)
            print(length(hvg))
            if(length(hvg) > 0){
                sobj <- subset(sobj, features=hvg)
            }
        } 
        con=runConos(sobj, batch)
        print(con)
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
        print(tool)
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







