#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 25.11.2022
Title: processes.nf

Description:
 - This script contains nextflow processes that are used to execute the scripts
 - The processes controls the activation of the desired conda environments, passes the input parameters and executes the scripts, and retunrs the output / writes the output to file
 - Workflow.nf script controls the execution of these scripts
"""

envDir = "$HOME/.conda/envs" // Requires miniconda3 in home dir


// ############################################################################### PREPARE DATA FOR INTEGRATION ALGORITHMS ##############################################################################################


//Process that returns list of methods to be benchmarked
process parse_methods {
    //Conda environment to be used
    conda "$envDir/PYenv"

    output: //Outputs the standard output of the executed process.
        stdout
    //Python script which, read the values of tools from config file, if the tool has value true, add it to list. After script is executed the list is returned
    script:
    """
        #!/usr/bin/env python
        tools_to_benchmark = []
        tools = {"scanvi":"$params.benchmark_scanvi", "scvi":"$params.benchmark_scvi", "scanorama":"$params.benchmark_scanorama", "combat":"$params.benchmark_combat", "mnn":"$params.benchmark_mnn", "desc":"$params.benchmark_desc", "saucie":"$params.benchmark_saucie", "bbknn":"$params.benchmark_bbknn", "harmony":"$params.benchmark_harmony", "scgen":"$params.benchmark_scgen", "trvae":"$params.benchmark_trvae", "conos":"$params.benchmark_conos", "fastmnn":"$params.benchmark_fastmnn", "seurat_cca":"$params.benchmark_seurat_cca","seurat_rpca":"$params.benchmark_seurat_rpca", "liger":"$params.benchmark_liger"}
        for tool,value in tools.items():
            if value == "true":
                tools_to_benchmark.append(tool)
        print(tools_to_benchmark)
    """
}

//Process that runs data processing, and returns the preprocessed adata files, the executed preprocessing notebooks for each method to be benchmarked, and .csv file containing paths to the preprocessed adata for each method to be benchmarked
process preprocessing{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Preprocessing/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.{h5ad, rds}" 
    //publishDir "Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.{h5ad}" 

    input: //Input variables
        path(input_data)
        path(benchmarking_settings_file)
        val(tools_to_benchmark)

    output: //Output variables
        path("*_*.h5ad") optional true
        path("*_*.rds") optional true
        path("*_*.RDS") optional true
        path("preprocessed_adata.csv", emit:preprocessed_adata_csv)
        

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/run_preprocessing.py \\
            $input_data \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings_file \\
            $tools_to_benchmark 
            > run_preprocessing_notebook.txt 2>&1
    """
}




// ########################################################################## INTEGRATION ALGORITHMS ########################################################################################################


process run_scvi{
    //Conda environment to be used
    conda "$envDir/PYenv"

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    input:
        val(algorithm)
        path(preprocessed_adata_csv)
        val(benchmarking_state)
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val(true, emit:benchmarking_state)
        //path("scvi_adata.h5ad")

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_scanvi{
    //Conda environment to be used
    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    conda "$envDir/PYenv"

    input:
        val(algorithm)
        path(preprocessed_adata_csv)
        val(benchmarking_state)
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val(true, emit:benchmarking_state)

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_harmony{
    //Conda environment to be used
    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    conda "$envDir/harmony_env"

    input:
        val(algorithm)
        path(preprocessed_adata_csv)
        val(benchmarking_state)
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val(true, emit:benchmarking_state)

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_saucie{
    //Conda environment to be used

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    conda "$envDir/saucie_env"

    input:
        val(algorithm)
        path(preprocessed_adata_csv)
        val(benchmarking_state)
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val(true, emit:benchmarking_state)

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_scanorama{
    //Conda environment to be used

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    conda "$envDir/scanorama_env"

    //Create the output dir for the process
    //publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state
        //path("scanorama_adata.h5ad")

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}


process run_combat{
    //Conda environment to be used

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    conda "$envDir/combat_env"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state
        //path("combat_adata.h5ad")


    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_liger{
    //Conda environment to be used

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    conda "$envDir/PYenv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}


process run_desc{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    //Conda environment to be used
    conda "$envDir/desc_env"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}


process run_bbknn{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    //Conda environment to be used
    conda "$envDir/bbknn_env"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_scgen{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_trvae{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_mnn{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true, pattern: "*.h5ad"

    //Conda environment to be used
    conda "$envDir/mnn_env"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)

    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/launch_integration_algorithms.py \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_seurat_cca{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true
    

    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)
    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    Rscript $projectDir/Scripts/launch_R_integration_algorithms.R \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $projectDir/Scripts/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_seurat_rpca{

    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true

    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)
    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    Rscript $projectDir/Scripts/launch_R_integration_algorithms.R \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $projectDir/Scripts/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_fastmnn{
    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true

    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)
    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    Rscript $projectDir/Scripts/launch_R_integration_algorithms.R \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $projectDir/Scripts/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}

process run_conos{
    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true

     //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    input:
        val algorithm
        path(preprocessed_adata_csv)
        val benchmarking_state
        path(benchmarking_settings)
    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    Rscript $projectDir/Scripts/launch_R_integration_algorithms.R \\
            $preprocessed_adata_csv \\
            $algorithm \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $projectDir/Scripts/ \\
            $benchmarking_settings
            > run_preprocessing_notebook.txt 2>&1
    """
}


// ######################################################################### INTEGRATION RESULT ################################################################################################

process convert_rds_output{
    publishDir "Output/Integration/Integrated_adata/", mode:"copy", overwrite:true

     //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    input:
        val algorithm
        val benchmarking_state
    
    output: //Outputs the state of the process
        val true, emit:benchmarking_state

    script: //Script to be executed
    """
    python $projectDir/Scripts/convert_rds_output.py \\
            $baseDir/Output/Integration/Integrated_adata/ \\
            $algorithm
            > run_preprocessing_notebook.txt 2>&1
    """
}

process parametrize_plotting_notebook{
    //Conda environment to be used
    conda "$envDir/PYenv"
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    //Create the output dir for the process
    //publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val tools
        val benchmarking_state
        path(benchmarking_settings_file)
    output:
        path("*.ipynb") 
        val(true, emit:process_state)


    script: //Script to be executed
    """
    python3 $projectDir/Scripts/plotting_notebook_parameters.py \\
            $baseDir/Output/Integration/Integrated_adata \\
            $benchmarking_settings_file \\
            $tools
            > run_preprocessing_notebook.txt 2>&1
    """
}


process run_results_notebook{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val process_state 

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/integration_results.ipynb --output integration_results.ipynb
    """
}