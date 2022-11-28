"""
Author: Erno Hänninen
Created: 2022-11-23
title main.nf

Description:
- connected to pipeline.config file
- controls which processes the program executes
- calls the workflow.nf file
"""


nextflow.enable.dsl=2


script_folder = "$baseDir/Scripts"

//Get processes from processes.nf file
//include {wf_parse_methods} from "$script_folder/workflow.nf"
include {wf_preprocessing} from "$script_folder/workflow.nf"
include {wf_integration_parameters} from "$script_folder/workflow.nf"
include {wf_integration_benchmarking} from "$script_folder/workflow.nf"



/*include {wf_scanvi} from "$script_folder/workflow.nf"
include {wf_scvi} from "$script_folder/workflow.nf"
include {wf_scgen} from "$script_folder/workflow.nf"
include {wf_scanorama} from "$script_folder/workflow.nf"
include {wf_bbknn} from "$script_folder/workflow.nf"
include {wf_combat} from "$script_folder/workflow.nf"
include {wf_desc} from "$script_folder/workflow.nf"
include {wf_harmony} from "$script_folder/workflow.nf"
include {wf_mnn} from "$script_folder/workflow.nf"
include {wf_saucie} from "$script_folder/workflow.nf"
include {wf_trvae} from "$script_folder/workflow.nf"
include {wf_trvaep} from "$script_folder/workflow.nf"*/



workflow {
    if(params.run_benchmarking){
        //workflow which parse which methods are benchmarked
        
        wf_preprocessing(params.input_data, params.benchmarking_settings)
        wf_integration_parameters(wf_preprocessing.out.preprocessed_adata_csv, params.benchmarking_settings)
        wf_integration_benchmarking(wf_integration_parameters.out.integration_parameters_state)
        /*if(!params.run_scanvi){
            wf_scanvi()
        }     
        if(!params.run_scvi){
            wf_scvi()
        }
        if(!params.run_scgen){
            wf_scgen()
        }
        if(!params.benchmark_scanorama){
            wf_scanorama()
        }
        if(!params.bbknn){
            wf_bbknn()
        }
        if(!params.combat){
            wf_combat()
        }
        if(!params.desc){
            wf_desc()
        }
        if(!params.harmony){
            wf_harmony()
        }
        if(!params.mnn){
            wf_mnn()
        }
        if(!params.saucie){
            wf_saucie()
        }
        if(!params.trvae){
            wf_trvae()
        }
        if(!params.trvaep){
            wf_trvaep()
        }*/

    }
}