#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 2022-11-23
Title: processes.nf

Description:
 - This script contains nextflow processes that are used to execute the scripts
 - The processes controls the activation of the desired conda environments, passes the input parameters and executes the scripts, and returns the output
"""

envDir = "$HOME/.conda/envs" // Requires miniconda3 in home dir


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

//Process that updates the parameters for notebooks to be executed in benchmerking.
/*process integration_parameters{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "Output/Integration/Integrated_adata", mode:"copy", overwrite:true, pattern: "*.h5ad"
    publishDir "Output/Integration/Metrics", mode:"copy", overwrite:true, pattern:".csv"


    //publishDir "Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.h5ad"
    input: //Input variables
        path(preprocessed_adata_csv)
        path(benchmarking_settings)

    output: //Output variables
        path("*.ipynb") 
        val true, emit:process_state

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/integration_notebook_parameters.py \\
            $preprocessed_adata_csv \\
            $benchmarking_settings
            > integration_notebook_parameters.txt 2>&1
    """
}

process run_scanvi{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val process_state
        val benchmarking_state

    
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_scanvi.ipynb --output run_scanvi.ipynb
    """
}

process run_scvi{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val process_state 
        val benchmarking_state

    
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_scvi.ipynb --output run_scvi.ipynb
    """
}

process run_scanorama{
    //Conda environment to be used
    conda "$envDir/scanorama_env"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    
    input:
        val process_state
        val benchmarking_state

    
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_scanorama.ipynb --output run_scanorama.ipynb
    """
}

process run_combat{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val process_state
        val benchmarking_state

        
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_combat.ipynb --output run_combat.ipynb
    """
}

process run_mnn{
    //Conda environment to be used
    conda "$envDir/mnn_env"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state

        
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_mnn.ipynb --output run_mnn.ipynb
    """
}

process run_desc{
    //Conda environment to be used
    conda "$envDir/desc_env"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state
    
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_desc.ipynb --output run_desc.ipynb
    """
}

process run_saucie{
    //Conda environment to be used
    conda "$envDir//saucie_env"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val process_state
        val benchmarking_state

    
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_saucie.ipynb --output run_saucie.ipynb
    """
}



process run_harmony{
    //Conda environment to be used
    conda "$envDir/harmony_env"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val process_state
        val benchmarking_state

    
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_harmony.ipynb --output run_harmony.ipynb
    """
}


process run_trvae{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state

        
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_trvae.ipynb --output run_trvae.ipynb
    """
}

process run_trvae_tl{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state

            
    output: //Outputs the state of the process
        val true, emit:process_state

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_trvae_tl.ipynb --output run_trvae_tl.ipynb
    """
}

process run_conos{
    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state


    output: //Output variables, the conos_state variable is used to trigger the quantify_conos process
        val true, emit:conos_state
        

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_conos.ipynb --output run_conos.ipynb
    """
}

process quantify_conos{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val conos_state
    
    output: //Outputs the state of the process
        val true, emit:process_state
        
    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/quantify_conos.ipynb --output quantify_conos.ipynb
    """
}

process run_fastmnn{
    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state

    output: //Output variables, the fastmnn_state variable is used to trigger the quantify_fastmnn process
        val true, emit:fastmnn_state
        

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_fastmnn.ipynb --output run_fastmnn.ipynb
    """
}

process quantify_fastmnn{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val fastmnn_state
    
    output: //Outputs the state of the process
        val true, emit:process_state
        
    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/quantify_fastmnn.ipynb --output quantify_fastmnn.ipynb
    """
}

process run_seurat_cca{
    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state

    output: //Output variables, the seurat_cca_state variable is used to trigger the quantify_seurat_cca process
        val true, emit:seurat_cca_state
        

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_seurat_cca.ipynb --output run_seurat_cca.ipynb
    """
}

process quantify_seurat_cca{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val seurat_cca_state
    
    output: //Outputs the state of the process
        val true, emit:process_state
        
    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/quantify_seurat_cca.ipynb --output quantify_seurat_cca.ipynb
    """
}


process run_seurat_rpca{
    //Conda environment to be used
    conda "$envDir/Renv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    input:
        val process_state
        val benchmarking_state

    output: //Output variables, the seurat_rpca_state variable is used to trigger the quantify_seurat_rpca process
        val true, emit:seurat_rpca_state
        

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/run_seurat_rpca.ipynb --output run_seurat_rpca.ipynb
    """
}

process quantify_seurat_rpca{
    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"

    input:
        val seurat_rpca_state
    
    output: //Outputs the state of the process
        val true, emit:process_state
        
    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/Output/Integration/Notebooks/quantify_seurat_rpca.ipynb --output quantify_seurat_rpca.ipynb
    """
}



process collect_integration_metrics{

    //Conda environment to be used
    conda "$envDir/PYenv"

    //Create the output dir for the process
    publishDir "Output/Integration/Metrics", mode:"copy", overwrite:true, pattern:"*.csv"

    input:
        val tools
        val integration_state

    output: //Output variables
        path("*.csv")
        
    script: //Script to be executed
    """
    python3 $projectDir/Scripts/collect_integration_metrics.py \\
            $tools 
            > integration_notebook_parameters.txt 2>&1
    """
}*/