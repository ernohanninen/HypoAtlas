#!/usr/bin/env nextflow
"""
Author: Erno Hänninen
Created: 25.11.2022
title main.nf

Description:
- calls the workflows in workflow.nf file
"""

nextflow.enable.dsl=2

script_folder = "$baseDir/Scripts"

//Get processes from processes.nf file
include {wf_preprocessing} from "$script_folder/workflow.nf"
include {wf_integration_benchmarking} from "$script_folder/workflow.nf"
include {wf_run_metrics} from "$script_folder/workflow.nf"

workflow {
    if(params.run_benchmarking){
        //workflow which parse which methods are benchmarked
        wf_preprocessing(params.input_data_path, params.benchmarking_settings)

        //Workflow which runs the integration algorithms
        wf_integration_benchmarking(wf_preprocessing.out.preprocessed_adata_csv, params.benchmarking_settings)

        //workflow which quantifies the integration algorithms
        wf_run_metrics(wf_preprocessing.out.tools_to_benchmark, wf_integration_benchmarking.out.benchmarking_state, params.benchmarking_settings)
    }
}