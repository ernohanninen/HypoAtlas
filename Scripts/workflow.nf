#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 2022-11-23
title workflow.nf

Description:
- This script works between main.nf and processes.nf scripts
"""

nextflow.enable.dsl=2

script_folder = "$baseDir/Scripts"

//Import the process modules from process.nf
include {parse_methods} from "$script_folder/processes.nf"
include {preprocessing} from "$script_folder/processes.nf"
include {integration_parameters} from "$script_folder/processes.nf"
include {run_scanvi} from "$script_folder/processes.nf"
include {run_scvi} from "$script_folder/processes.nf"
include {run_scanorama} from "$script_folder/processes.nf"
include {run_combat} from "$script_folder/processes.nf"
include {run_mnn} from "$script_folder/processes.nf"
include {run_desc} from "$script_folder/processes.nf"
include {run_saucie} from "$script_folder/processes.nf"
include {run_bbknn} from "$script_folder/processes.nf"
include {run_harmony} from "$script_folder/processes.nf"
include {run_scgen} from "$script_folder/processes.nf"
include {run_trvae} from "$script_folder/processes.nf"
include {run_trvae_tl} from "$script_folder/processes.nf"
include {run_conos} from "$script_folder/processes.nf"
include {quantify_conos} from "$script_folder/processes.nf"
include {run_fastmnn} from "$script_folder/processes.nf"
include {quantify_fastmnn} from "$script_folder/processes.nf"
include {run_liger} from "$script_folder/processes.nf"
include {quantify_liger} from "$script_folder/processes.nf"
include {run_seurat_cca} from "$script_folder/processes.nf"
include {quantify_seurat_cca} from "$script_folder/processes.nf"
include {run_seurat_rpca} from "$script_folder/processes.nf"
include {quantify_seurat_rpca} from "$script_folder/processes.nf"
include {collect_integration_metrics} from "$script_folder/processes.nf"



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

//Workflow which calls the process that is used to adjust the parameters of integration notebooks
workflow wf_integration_parameters{
    take: //Workflow input
        preprocessed_adata_csv
        benchmarking_settings
    main: //Calling processes
        integration_parameters(preprocessed_adata_csv, benchmarking_settings)
    emit: //Output
        integration_parameters_state = integration_parameters.out.process_state
}

workflow wf_integration_benchmarking{
    take: //Input
        tools
        integration_parameters_state
    main: //Calling the processes

        //Condititional workflow, the params. -parameters are readed from pipeline.config file
        //The state variable makes the wf_collect_metrics workflow dependent of wf_integration_benchmarking workflow, therefore preventing the execution of wf_collect_metrics before this workflow is ready
        if(params.benchmark_scanvi){ state = run_scanvi(integration_parameters_state) }
        if(params.benchmark_scvi){ state = run_scvi(integration_parameters_state) }
        if(params.benchmark_scanorama) { state = run_scanorama(integration_parameters_state) }
        if(params.benchmark_combat) { state = run_combat(integration_parameters_state) }
        if(params.benchmark_mnn) { state = run_mnn(integration_parameters_state) }
        if(params.benchmark_desc) { state = run_desc(integration_parameters_state) }
        if(params.benchmark_saucie) { state = run_saucie(integration_parameters_state) }
        if(params.benchmark_bbknn) { state =  run_bbknn(integration_parameters_state) }
        if(params.benchmark_harmony) { state = run_harmony(integration_parameters_state) }
        if(params.benchmark_scgen) { state = run_scgen(integration_parameters_state) }
        if(params.benchmark_trvae) { state = run_trvae(integration_parameters_state) }
        if(params.benchmark_trvae_tl) { state = run_trvae_tl(integration_parameters_state) }
        //The rest of the integration methods are written in R, therefore there is an external process running the quantification of the method (written in python)
        if(params.benchmark_conos) {
            run_conos(integration_parameters_state) 
            state = quantify_conos(run_conos.out.conos_state)  //The run_conos.out.conos_state state variable takes care that the wuantification process is not executed before the run_conos is ready
        }
        if(params.benchmark_fastmnn) {
            run_fastmnn(integration_parameters_state)
            state = quantify_fastmnn(run_fastmnn.out.fastmnn_state)   
        }
        if(params.benchmark_liger) {
            run_liger(integration_parameters_state)
            state = quantify_liger(run_liger.out.liger_state)
        }

        if(params.benchmark_seurat_cca) {
            run_seurat_cca(integration_parameters_state)
            state = quantify_seurat_cca(run_seurat_cca.out.seurat_cca_state)
        }

        if(params.benchmark_seurat_rpca) {
            run_seurat_rpca(integration_parameters_state)
            state = quantify_seurat_rpca(run_seurat_rpca.out.seurat_rpca_state)
        }

    emit: //Outputting the state variable that is updated in every if statement
        benchmarking_state = state

}

//Workflow which calls process that collects the integration metrics of different algorithms to one file
workflow wf_collect_metrics{
    take: //Input
        tools
        benchmarking_state
    main: //Calling the process
        collect_integration_metrics(tools, benchmarking_state)

}

// ######################################################################################################################################################
