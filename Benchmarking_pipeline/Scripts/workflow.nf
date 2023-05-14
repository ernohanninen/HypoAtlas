#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 25.11.2023
title workflow.nf

Description:
- This script works between main.nf and processes.nf scripts
- Contains the workflows that calls the processes
- PAsses arguments to the processes

Workflow procedure:
- Preparing (find out the methods to be benchmarked, run preprocessing)
- Workflow for integration benchmarking (Involves running all the methods to be benchmarked)
- Workflow for integration results (parameterize notebook and execute it)

"""

nextflow.enable.dsl=2

script_folder = "$baseDir/Scripts"

//Import the process modules from process.nf
include {parse_methods} from "$script_folder/processes.nf"
include {preprocessing} from "$script_folder/processes.nf"
include {run_scvi} from "$script_folder/processes.nf"
include {run_scanvi} from "$script_folder/processes.nf"
include {run_scanorama} from "$script_folder/processes.nf"
include {run_combat} from "$script_folder/processes.nf"
include {run_harmony} from "$script_folder/processes.nf"
include {run_liger} from "$script_folder/processes.nf"
include {run_desc} from "$script_folder/processes.nf"
include {run_scgen} from "$script_folder/processes.nf"
include {run_bbknn} from "$script_folder/processes.nf"
include {run_seurat_cca} from "$script_folder/processes.nf"
include {run_seurat_rpca} from "$script_folder/processes.nf"
include {convert_rds_output} from "$script_folder/processes.nf"
include {run_fastmnn} from "$script_folder/processes.nf"
include {run_conos} from "$script_folder/processes.nf"
include {run_trvae} from "$script_folder/processes.nf"
include {run_mnn} from "$script_folder/processes.nf"
include {run_saucie} from "$script_folder/processes.nf"
include {parametrize_plotting_notebook} from "$script_folder/processes.nf"
include {run_results_notebook} from "$script_folder/processes.nf"

//######################################################### INTEGRATION BENCHMARKING WORKFLOW ###################################################################

//Workflwo which first calls the parse_methods() process to get the methods to be benchmarked, then calls then preprocessing() process that runs the preprocessing 
workflow wf_preprocessing{

    take: //Input of the workflow
        input_data //change this to h5ad object, now it is the file with all the data
        benchmarking_settings
    
    main: //Calling the processes from processes.nf
        tools_to_benchmark = parse_methods()
        preprocessing(input_data, benchmarking_settings, tools_to_benchmark)
    
    emit: //Output of the process
        preprocessed_adata_csv = preprocessing.out.preprocessed_adata_csv  
        tools_to_benchmark = tools_to_benchmark     
}

workflow wf_integration_benchmarking{
    take: //Input
        preprocessed_adata_csv
        benchmarking_settings
        
    main: //Calling the processes
        benchmarking_state = false
        //Condititional workflow
        //The state variable makes the wf_collect_metrics workflow dependent of wf_integration_benchmarking workflow, therefore preventing the execution of wf_collect_metrics before this workflow is ready
        if(params.benchmark_scvi) { benchmarking_state = run_scvi("scvi", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_scanvi) { benchmarking_state = run_scanvi("scanvi", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_scanorama) { benchmarking_state = run_scanorama("scanorama", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_combat) { benchmarking_state = run_combat("combat", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_harmony) { benchmarking_state = run_harmony("harmony", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_liger) { benchmarking_state = run_liger("liger", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_desc) { benchmarking_state = run_desc("desc", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_bbknn) { benchmarking_state = run_bbknn("bbknn", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_scgen) { benchmarking_state = run_scgen("scgen", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_trvae) { benchmarking_state = run_trvae("trvae", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_mnn) { benchmarking_state = run_mnn("mnn", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }
        if(params.benchmark_saucie) { benchmarking_state = run_saucie("saucie", preprocessed_adata_csv, benchmarking_state, benchmarking_settings) }

        if(params.benchmark_conos) {
            benchmarking_state = run_conos("conos", preprocessed_adata_csv, benchmarking_state, benchmarking_settings)

        }
        if(params.benchmark_fastmnn) {
            benchmarking_state = run_fastmnn("fastmnn", preprocessed_adata_csv, benchmarking_state, benchmarking_settings)
        }

        if(params.benchmark_seurat_cca) {
            benchmarking_state = run_seurat_cca("seurat_cca", preprocessed_adata_csv, benchmarking_state, benchmarking_settings)
        }

        if(params.benchmark_seurat_rpca) {
            benchmarking_state = run_seurat_rpca("seurat_rpca", preprocessed_adata_csv, benchmarking_state, benchmarking_settings)

        }
    emit: //Outputting the state variable that is updated in every if statement
        benchmarking_state

}

// Workflow for integration results
workflow wf_run_metrics{
    take: //Input
        tools
        benchmarking_state
        benchmarking_settings       
        
    main: //Calling the processes
        parametrize_plotting_notebook(tools, benchmarking_state, benchmarking_settings)
        run_results_notebook(parametrize_plotting_notebook.out.process_state)
}
