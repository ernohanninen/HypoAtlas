#!/usr/bin/env python
"""
Author: Erno Hänninen
Created: 2022-11-24
Title: run_preprocessing_notebook.py

Description:
    - Read notebook parameters from benchmarking_settings.csv -> update notebook -> run preprocessing notebook 

Procedure:
    - Loop over list, that contains the tools to be benchmarked
    - During every iteration:
        1. Extract the tool specific parameters from benchmarking_settings.csv file
        2. Update the parameters to the run_preprocessing.ipynb notebook
        3. Run the preprocessing notebook, and store the executed notebook to output folder
            - The run_preprocessing.ipynb notebook performs data processing and stores the processed adata object to output folder
            
List of functions:
    - extract_preprocessing_settings, update_notebook_parameters, run_preprocessing

List of non-standard modules:
    - pandas, nbclient, nbformat, nbparameterise
"""


import scanpy as sc
import subprocess #Allows running bash
import os, sys, csv, re, shutil
import pandas as pd


#Function which extracts the preprocessing settings for each tool
def extract_preprocessing_settings(tool, settings_df):
    #Get the settings of specific tool and read them to variable
    method, batch, num_hvg, r_output, seurat_output, scale, label_key = settings_df.loc[settings_df["integration_method"] == tool].values[0]
    
    # Reformat the input values from csv file
    if num_hvg == "None" or num_hvg == None: #if num_hvg is none convert it to 0. In the preprocessing script if num_hvg is < 500, the hvg is not computed -> meaning that the all genes in dataset are used in integration
        num_hvg = 0 
    num_hvg = int(num_hvg) 
    
    #Function that converts the string true and false values to boolean
    def str_to_bool(s): 
        if s == 'True':
            return True
        elif s == 'False':
            return False
        else:
            return s
    r_output = str_to_bool(r_output)
    seurat_output = str_to_bool(seurat_output)
    scale = str_to_bool(scale)
    
    print("SETTINGS TO BE USED:")
    print(method, batch, num_hvg, r_output, seurat_output, scale, label_key)
    return  method, batch, num_hvg, r_output, seurat_output, scale, label_key #return settings

#Function which updates notebook parameters
"""def update_notebook_parameters(input_data_path, integration_method, batch, num_hvg, r_output, seurat_output, scale):
    #Open the preprocessing notebook
    with open("../../../BenchmarkStudy/Notebooks_to_run/run_preprocessing.ipynb") as f:
        nb = nbformat.read(f, as_version=4) #read the notebook to variable

    # Get a list of Parameter objects
    orig_parameters = extract_parameters(nb)
    

    # Update one or more parameters
    params = parameter_values(orig_parameters, method=integration_method,  file=os.getcwd() + "/" + input_data_path,  out = os.getcwd(), hvg=num_hvg, rout=r_output, batch=batch, seurat = seurat_output, scale = scale)

    # Make a notebook object with these definitions
    new_nb = replace_definitions(nb, params)

    #Write the edited notebook to file
    nbformat.write(new_nb, fp="../../../BenchmarkStudy/Notebooks_to_run/run_preprocessing.ipynb")
"""  
    
def run_preprocessing(input_data_path, method, batch, num_hvg, scale, r_output, seurat_output):
    adata = sc.read(input_data_path)
    hvgs = adata.var.index
    
    print("input_data_path: ", input_data_path)
    print("method: ", method)
    print("batch: ", batch)
    print("hvg: ", num_hvg)
    print("scale: ", scale)
    print("r_output : ", r_output)
    print("seurat_output: ", seurat_output )

    # remove HVG if already precomputed
    if 'highly_variable' in adata.var:
        del adata.var['highly_variable']

    if num_hvg > 500:
        print("Computing HVGs ...")
        if seurat_output == True:
            hvgs = scib.preprocessing.hvg_batch(
                adata,
                batch_key=batch,
                target_genes=num_hvg,
                adataOut=False
            )
        else:
            adata = scib.preprocessing.hvg_batch(
                adata,
                batch_key=batch,
                target_genes=num_hvg,
                adataOut=True
            )
    
    
    if scale == True:
        print("Scaling data ...")
        adata = scib.preprocessing.scale_batch(adata, batch)
   
    if r_output == True:
        print("Save as RDS")
        scib.preprocessing.saveSeurat(adata, os.getcwd() + "/"+ method + "_seurat.rds", batch, hvgs)

    else:
        print("Save as HDF5")
        adata.write(os.getcwd() + "/"+ method + "_adata.h5ad")
    
    #Depending of the output type return the correct output path
    if r_output == True:
        return os.getcwd() + "/" + method + "_seurat.rds"
    return os.getcwd() + "/" + method + "_adata.h5ad"
    
    
#This if statement is executed when script is called
if __name__ == "__main__":

    import scib
    # packages which allows edit notebooks from python
    import nbclient
    import nbformat
    from nbparameterise import (
        extract_parameters, replace_definitions, parameter_values
    )

    #Read input arguments
    input_data_path = sys.argv[1]
    output_data_path = sys.argv[2]
    benchmarking_settings_file = sys.argv[3]
    tools_to_benchmark = sys.argv[4:int(len(sys.argv))]
    
    #Read the benchmarking_settings.csv file to pandas df
    settings_df = pd.read_csv(benchmarking_settings_file)
    print(settings_df)
    
    adata_dict = {} #store the processed adata file paths to dict
    
    if os.path.exists(output_data_path):
        shutil.rmtree(output_data_path)
    os.makedirs(output_data_path)
    
    #loop over the tools, during every iteration call the extract_preprocessing_settings and run_preprocessing functions
    for tool in tools_to_benchmark: 
        tool = re.sub(r'[\W*]', '', tool) #As nextflow messes the python list created in parse_methods process, remove the unwanted '[', ',' and ']' -characters from the list items 
        method, batch, num_hvg, r_output, seurat_output, scale, label_key = extract_preprocessing_settings(tool, settings_df)  #Extract settings from the input df
        file_path = run_preprocessing(input_data_path, method, batch, num_hvg, scale, r_output, seurat_output) #Update parameters to notebook
        
        #file_path = run_preprocessing(method, r_output) #Runs preprocessing 
        adata_dict.update({tool:file_path}) #Store the preprocessed adata file path to dict 
    
    adata = sc.read(input_data_path)
    adata = scib.preprocessing.hvg_batch(adata,batch_key=batch,target_genes=num_hvg,adataOut=True)
    sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
    
    
    adata.obsm["unintegrated"] = adata.obsm["X_pca"]
    
    print("ADATA UNINTEGRATED: ", adata)
    print("adata.X: ", adata.X)
    print("adata.layers[counts]: ", adata.layers["counts"])

    adata.write(output_data_path + "unintegrated.h5ad")
        
    
    #Store the file paths to csv file
    with open("preprocessed_adata.csv", "w") as csvfile:
        fieldnames= ["integration_method", "adata_path"] #Initialize the headers
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames) #Create writer object
        writer.writeheader()
        for key, value in adata_dict.items(): #Loop over the dict and store the key and value to file
            writer.writerow({"integration_method":key, "adata_path":value}) 
        writer.writerow({"integration_method":"unintegrated", "adata_path":os.getcwd() + "/unintegrated.h5ad"})    
           