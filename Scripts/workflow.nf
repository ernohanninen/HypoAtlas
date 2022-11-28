//This script works between main.nf and processes.nf scripts

nextflow.enable.dsl=2

script_folder = "$baseDir/Scripts"
include {parse_methods} from "$script_folder/processes.nf"
include {preprocessing} from "$script_folder/processes.nf"
include {integration_parameters} from "$script_folder/processes.nf"
include {run_scanvi} from "$script_folder/processes.nf"
include {run_scvi} from "$script_folder/processes.nf"
include {run_scanorama} from "$script_folder/processes.nf"

/*include {run_scgen} from "$script_folder/workflow.nf"
include {run_scanorama} from "$script_folder/workflow.nf"
include {run_bbknn} from "$script_folder/workflow.nf"
include {run_combat} from "$script_folder/workflow.nf"
include {run_desc} from "$script_folder/workflow.nf"
include {run_harmony} from "$script_folder/workflow.nf"
include {run_mnn} from "$script_folder/workflow.nf"
include {run_saucie} from "$script_folder/workflow.nf"
include {run_trvae} from "$script_folder/workflow.nf"
include {run_trvaep} from "$script_folder/workflow.nf*/


/*include {run_scanvi} from "$script_folder/processes.nf"
include {run_scvi} from "$script_folder/processes.nf"
include {run_scgen} from "$script_folder/processes.nf"*/


workflow wf_preprocessing{
    take:
        input_data //change this to h5ad object, no it is the file with all the data
        benchmarking_settings
    
    main:
        tools_to_benchmark = parse_methods()
        preprocessing(input_data, benchmarking_settings, tools_to_benchmark)
    
    emit:
        preprocessed_adata_csv = preprocessing.out.preprocessed_adata_csv       
}


workflow wf_integration_parameters{
    take:
        preprocessed_adata_csv
        benchmarking_settings
    main:
        integration_parameters(preprocessed_adata_csv, benchmarking_settings)
    emit:
        integration_parameters_state = integration_parameters.out.process_state


}

workflow wf_integration_benchmarking{
    take:
        integration_parameters_state

    main:
        if(params.benchmark_scanvi){ run_scanvi(integration_parameters_state) }
        if(params.benchmark_scvi){ run_scvi(integration_parameters_state) }
        if(params.benchmark_scanorama) { run_scanorama(integration_parameters_state )}
        /*if(!params.run_scgen){
            run_scgen()
        }
        if(!params.benchmark_scanorama){
            run_scanorama()
        }
        if(!params.bbknn){
            run_bbknn()
        }
        if(!params.combat){
            run_combat()
        }
        if(!params.desc){
            run_desc()
        }
        if(!params.harmony){
            run_harmony()
        }
        if(!params.mnn){
            run_mnn()
        }
        if(!params.saucie){
            run_saucie()
        }
        if(!params.trvae){
            run_trvae()
        }
        if(!params.trvaep){
            run_trvaep()
        }*/
}



