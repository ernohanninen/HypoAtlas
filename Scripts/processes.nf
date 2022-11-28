#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 2022-11-23
Title: processes.nf

Description:
 - This script calls the scripts
"""


//Process that returns list of methods to be benchmarked
process parse_methods {
    //Conda environment to be used
    conda "/home/ernohanninen/miniconda3/envs/PYenv"

    output: //Outputs the standard output of the executed process.
        stdout
    //Python script which, read the values of tools from config file, if the tool has value true, add it to list. After script is executed the list is returned
    script:
    """
        #!/usr/bin/env python
        tools_to_benchmark = []
        tools = {"scanvi":"$params.benchmark_scanvi", "scvi":"$params.benchmark_scvi", "scgen":"$params.benchmark_scgen", "scanorama":"$params.benchmark_scanorama"}
        for tool,value in tools.items():
            if value == "true":
                tools_to_benchmark.append(tool)
        print(tools_to_benchmark)
    """
}

//Process that runs data processing, and returns the preprocessed adata files, the executed preprocessing notebooks for each method to be benchmarked, and .csv file containing paths to the preprocessed adata for each method to be benchmarked
process preprocessing{
    //Conda environment to be used
    conda "/home/ernohanninen/miniconda3/envs/PYenv"

    //Create the output dir for the process
    publishDir "BenchmarkStudy/Output/Preprocessing/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    publishDir "BenchmarkStudy/Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.{h5ad, rds}" 
    //publishDir "BenchmarkStudy/Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.{h5ad}" 

    input: //Input variables
        path(input_data)
        path(benchmarking_settings_file)
        val(tools_to_benchmark)

    output: //Output variables
        path("*_*.h5ad") optional true
        path("*_*.rds") optional true
        path("*.ipynb")
        path("preprocessed_adata.csv", emit:preprocessed_adata_csv)

    script: //Script to be executed
    """
    python3 $projectDir/Scripts/run_preprocessing_notebook.py \\
            $input_data \\
            $benchmarking_settings_file \\
            $tools_to_benchmark 
            > run_preprocessing_notebook.txt 2>&1
    """
}

//Process that updates the parameters for notebooks to be executed in benchmerking.
process integration_parameters{
    //Conda environment to be used
    conda "/home/ernohanninen/miniconda3/envs/PYenv"

    //Create the output dir for the process
    publishDir "BenchmarkStudy/Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    //publishDir "BenchmarkStudy/Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.h5ad"

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
    conda "/home/ernohanninen/miniconda3/envs/PYenv"
    input:
        val process_state
    //Create the output dir for the process
    //publishDir "BenchmarkStudy/Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    //publishDir "BenchmarkStudy/Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.h5ad"
    output: //Output variables
        path("*.ipynb")

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/BenchmarkStudy/Notebooks_to_run/run_scanvi.ipynb --output scanvi_benchmarking.ipynb
    """
}

process run_scvi{
    //Conda environment to be used
    conda "/home/ernohanninen/miniconda3/envs/PYenv"

    //Create the output dir for the process
    //publishDir "BenchmarkStudy/Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    //publishDir "BenchmarkStudy/Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.h5ad"
    input:
        val process_state

    output: //Output variables
        path("*.ipynb")

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/BenchmarkStudy/Notebooks_to_run/run_scvi.ipynb --output scvi_benchmarking.ipynb
    """
}

process run_scanorama{
    //Conda environment to be used
    conda "/home/ernohanninen/miniconda3/envs/scanorama_env"

    //Create the output dir for the process
    //publishDir "BenchmarkStudy/Output/Integration/Notebooks", mode:"copy", overwrite:true, pattern: "*.ipynb"
    //publishDir "BenchmarkStudy/Output/Preprocessing/AnnData", mode:"copy", overwrite:true, pattern: "*.h5ad"
    input:
        val process_state

    output: //Output variables
        path("*.ipynb")

    script: //Script to be executed
    """
    jupyter nbconvert --to notebook --execute $projectDir/BenchmarkStudy/Notebooks_to_run/run_scanorama.ipynb --output scanorama_benchmarking.ipynb
    """
}