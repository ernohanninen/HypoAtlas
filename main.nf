#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 2022-11-23
title main.nf

Description:
- calls the workflow.nf file
"""

nextflow.enable.dsl=2

script_folder = "$baseDir/Scripts"

//Get processes from processes.nf file
include {wf_preprocessing} from "$script_folder/workflow.nf"
include {wf_integration_parameters} from "$script_folder/workflow.nf"
include {wf_integration_benchmarking} from "$script_folder/workflow.nf"
include {wf_collect_metrics} from "$script_folder/workflow.nf"

workflow {
    if(params.run_benchmarking){
        //workflow which parse which methods are benchmarked
        wf_preprocessing(params.input_data_path, params.benchmarking_settings)
        //Workflow which adjusts the integration parameters for different notebooks 
        wf_integration_parameters(wf_preprocessing.out.preprocessed_adata_csv, params.benchmarking_settings)
        //Workflow which runs the integration algorithms
        wf_integration_benchmarking(wf_preprocessing.out.tools_to_benchmark,wf_integration_parameters.out.integration_parameters_state)
        //Workflow which collects the integration metrics of different algorithms to one csv file
        wf_collect_metrics(wf_preprocessing.out.tools_to_benchmark, wf_integration_benchmarking.out.benchmarking_state)
    }
}