/*
Title: data_processing_wf.nf
Date: 15.12.2022
Author: Erno Hänninen
Description:
- Loads data from given url (if available) and processes the datasets required for this study

Workflow description:
- Read raw and processed herb data from given adress, call script which processes raw Herb data using processed data as reference
- Read raw zhou count matrices from GEO, the processed Zhou is read from local resources (not publicly available), and call which processes raw Zhou data using processed data as reference
- Call script which merges the zhou and herb data
*/


dataDir = "$baseDir/Data"
envDir = "$HOME/.conda/envs" // Requires miniconda3 in home dir
scriptDir = "$baseDir/Scripts"

// Process that converts the seurat object to h5ad file
process convertHerbData {
    // Load integrated and raw timepoint data (authors have calculated qc metric)
    // Note the Fig1_EmbryonicHypo_FINAL.RData is processed by the authors, but as it needs some processing it is loaded to Raw-folder
    beforeScript "wget -N -P $dataDir/RawHerb https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_AllCells.rda && wget -N -P $dataDir/RawHerb https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/Fig1_EmbryonicHypo_FINAL.RData"

    conda "$envDir/RConvertEnv" // Activate conda env

    input: // input is file path
        val raw_data
        val integrated_data

    output: // output files 
        path("MergedHerb.h5ad", emit: merged_raw_data)
        path("IntegratedHerb.h5ad", emit: integrated_data)

    script: // Script to execute
    """
        Rscript $scriptDir/convert_to_h5ad.R \\
        $dataDir/RawHerb/$raw_data \\
        $dataDir/RawHerb/$integrated_data 
        > convert_herb_data.txt 2>&1
    """
}

// Process that annotates the raw data based on the integrated object
process processHerbData {

    conda "$envDir/PYenv" 

    publishDir "Output", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "Data/Processed", mode:"copy", overwrite:true, pattern: "*.h5ad"

    input: 
        path(merged_raw_herb_h5ad)
        path(integrated_herb_h5ad)
        
    output: 
        path("Processed_herb_adata.h5ad", emit: processed_herb)
        path("*.ipynb")
        val(true, emit:herbState)

    script: //Script to execute
    // Parameterizes the notebook
    """
        papermill -p merged_raw_data $merged_raw_herb_h5ad -p integrated_data $integrated_herb_h5ad $scriptDir/AnnotateHerbData.ipynb AnnotateHerb_notebook.ipynb
    """

}

// Process that loads and annotates the zhou data
// Requires that an annotated reference is stored in the $dataDir/RawZhou
process processZhouData {

    // Loads the raw reads from GEO
    beforeScript "wget -N -P $dataDir/RawZhou ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169109/suppl/GSE169109_barcodes.tsv.gz && wget -N -P $dataDir/RawZhou https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169109/suppl/GSE169109_features.tsv.gz && wget -N -P $dataDir/RawZhou https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169109/suppl/GSE169109_matrix.mtx.gz"

    conda "$envDir/PYenv" // Activate conda env

    publishDir "Output", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "Data/Processed", mode:"copy", overwrite:true, pattern: "*.h5ad"

    input:
        val(herbState)

    output: 
        path("Processed_zhou_adata.h5ad", emit: processed_zhou)
        path("*.ipynb")

    script: 
    // Parameterizes  and executes the notebook

    """
        papermill -p annotated_data_path $dataDir/RawZhou/zhou_annotations.RDS -p raw_read_path $dataDir/RawZhou $scriptDir/AnnotateZhouData.ipynb AnnotateZhou_notebook.ipynb
    """
}

// Process that concatanates the annotated herb and zhou data
process mergeData {

    conda "$envDir/PYenv" // Activate conda env

    publishDir "Output", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "Data/Processed", mode:"copy", overwrite:true, pattern: "*.h5ad"

    input: 
        path(processedHerb)
        path(processedZhou)
        
    output: 
        path("merged_zhou_herb.h5ad")
        path("*.ipynb")

    script: 
    // Parameterizes  and executes the notebook

    """
        papermill -p herb_path $processedHerb -p zhou_path $processedZhou $scriptDir/Concat_Zhou_Herb.ipynb ConcatData_notebook.ipynb
    """

}


// Workflow controls the process execution
workflow { 
    convertHerbData("HypoSamples_AllCells.rda", "Fig1_EmbryonicHypo_FINAL.RData")
    processHerbData(convertHerbData.out.merged_raw_data, convertHerbData.out.integrated_data)
    processZhouData(processHerbData.out.herbState)
    mergeData(processHerbData.out.processed_herb, processZhouData.out.processed_zhou)

}