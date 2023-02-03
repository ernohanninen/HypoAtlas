



seurat_cca = function(data, batch, hvg=2000) {

	  batch_list = SplitObject(data, split.by = batch)

	  anchors = FindIntegrationAnchors(
            object.list = batch_list,
            anchor.features = hvg,
            scale = T,
            l2.norm = T,
            dims = 1:30,
            k.anchor = 5,
            k.filter = 200,
            k.score = 30,
            max.features = 200,
            eps = 0)

        integrated = IntegrateData(
            anchorset = anchors,
            new.assay.name = "integrated",
            features = NULL,
            features.to.integrate = NULL,
            dims = 1:30,
            k.weight = 100,
            weight.reduction = NULL,
            sd.weight = 1,
            sample.tree = NULL,
            preserve.order = F,
            eps = 0,
            verbose = T)
	  return(integrated)
}

# This function is from https://github.com/theislab/scib-pipeline/blob/0c7be53b1000864fcd31a7b7594f9a5071204233/scripts/integration/integration.R
#runSeuratRPCA = function(data, batch, hvg=2000) {
seuratRPCA = function(data, batch, hvg=2000) {
     
    batch_list = SplitObject(data, split.by = batch)
        #features <- SelectIntegrationFeatures(batch_list)
        batch_list <- lapply(X = batch_list, FUN = function(x) {
        x  <- ScaleData(x, features = hvg)
        x <- RunPCA(x, features = hvg)
        return(x)
    })

    anchors = FindIntegrationAnchors(
        object.list = batch_list,
        anchor.features = hvg,
        scale = T,
        l2.norm = T,
        dims = 1:30,
        k.anchor = 5,
        k.filter = 200,
        k.score = 30,
        reduction = "rpca",
        max.features = 200,
        eps = 0)
        integrated = IntegrateData(
        anchorset = anchors,
        new.assay.name = "integrated",
        features = NULL,
        features.to.integrate = NULL,
        dims = 1:30,
        k.weight = 100,
        weight.reduction = NULL,
        sd.weight = 1,
        sample.tree = NULL,
        preserve.order = F,
        eps = 0,
        verbose = T)
    return(integrated)
}

#This function is from https://github.com/theislab/scib-pipeline/blob/0c7be53b1000864fcd31a7b7594f9a5071204233/scripts/integration/integration.R
preP <- function(so, vars.to.regress=NULL, verbose=TRUE, n.pcs=100) {
  if (verbose) {
    message("Running Seurat v3 workflow")
  }
  so <- Seurat::FindVariableFeatures(object = so, verbose = verbose)
  so <- Seurat::ScaleData(object = so, verbose = verbose)
  so <- Seurat::RunPCA(object = so, npcs = n.pcs, verbose = verbose)
  return(so)
}


#This function is from https://github.com/theislab/scib-pipeline/blob/0c7be53b1000864fcd31a7b7594f9a5071204233/scripts/integration/integration.R
runConos = function(sobj, batch) {

    batch_list <- SplitObject(sobj, split.by=batch)
    pp <- lapply(batch_list, preP)
    con <- Conos$new(pp)
    con$buildGraph(space="genes")
    con$findCommunities()
    con$embedGraph(method="UMAP")
    
    return(con)
}


# This function is from https://github.com/theislab/scib-pipeline/blob/0c7be53b1000864fcd31a7b7594f9a5071204233/scripts/integration/integration.R
run_fastMNN = function(sobj, batch) {
    

    if (is.null(sobj@assays$RNA)) {
        # Seurat v4
        expr <- GetAssayData(sobj, slot = "data")
    } else {
        # Seurat v3
        expr <- sobj@assays$RNA@data
    }

    sce <- fastMNN(expr, batch = sobj@meta.data[[batch]])
    corrected_data <- assay(sce, "reconstructed")

    if (is.null(sobj@assays$RNA)) {
        # Seurat v4
        sobj <- SetAssayData(sobj, slot = "data", new.data = as.matrix(corrected_data))
        sobj@reductions['fastmnn'] <- CreateDimReducObject(reducedDim(sce, "corrected"), key = 'fastmnn_')
    } else {
        # Seurat v3
        sobj@assays$RNA <- CreateAssayObject(corrected_data)
        sobj[['fastmnn']] <- CreateDimReducObject(reducedDim(sce, "corrected"), key = 'fastmnn_')
    }
    print("INTEGRATION READY")

    return(sobj)
}
