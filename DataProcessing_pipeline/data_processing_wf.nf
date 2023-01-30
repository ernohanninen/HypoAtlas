
/*
Title: data_processing_wf.nf
Date: 2023-01-24
Author: Erno Hänninen
Description:
*/


dataDir = "$baseDir/Data"
envDir = "$HOME/.conda/envs" // Requires miniconda3 in home dir
scriptDir = "$baseDir/Scripts"


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

    script: // Script that converts the seurat object to h5ad file
    """
        Rscript $scriptDir/convert_to_h5ad.R \\
        $dataDir/RawHerb/$raw_data \\
        $dataDir/RawHerb/$integrated_data 
        > convert_herb_data.txt 2>&1
    """
}

process processHerbData {

    conda "$envDir/PYenv" // Activate conda env

    publishDir "Output", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "Data/Processed", mode:"copy", overwrite:true, pattern: "*.h5ad"


    input: 
        path(merged_raw_herb_h5ad)
        path(integrated_herb_h5ad)
        

    output: // output files 
        path("Processed_herb_adata.h5ad", emit: processed_herb)
        path("*.ipynb")



    script: // Script that converts the seurat object to h5ad file
    """
        papermill -p merged_raw_data $merged_raw_herb_h5ad -p integrated_data $integrated_herb_h5ad $scriptDir/AnnotateHerbData.ipynb AnnotateHerb_notebook.ipynb
    """

}

/*process loadZhouData {

    beforeScript "wget -N -P $dataDir/RawHerb https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_AllCells.rda && wget -N -P $dataDir/RawHerb https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/Fig1_EmbryonicHypo_FINAL.RData"


}

process processZhouData {

}

process mergeData {

}*/

workflow { 
    convertHerbData("HypoSamples_AllCells.rda", "Fig1_EmbryonicHypo_FINAL.RData")
    processHerbData(convertHerbData.out.merged_raw_data, convertHerbData.out.integrated_data)

    //mergeData(processHerbData.out.process_herb)

}