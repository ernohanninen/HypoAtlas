#!/usr/bin/env python
"""
Author: Erno Hänninen
Created: 2023-01-25
title run_integration_notebooks.py

Description:
- 
"""
import os
import scanpy as sc
import sys, csv
#import subprocess #Allows running bash
import pandas as pd
import re
#packages which allows edit notebooks from python
import nbclient, nbformat
from nbparameterise import (extract_parameters, replace_definitions, parameter_values)
#sys.path.insert(script_dir) #Adding a path to be able to import a function from run_preprocessing.py
from run_preprocessing import extract_preprocessing_settings
#import papermill as pm

def update_notebook_parameters(adata_unintegrated, integrated_adata, tools,data_dir, batch, label):
    #Open the preprocessing notebook
    with open("../../../Scripts/integration_results.ipynb") as f:
        nb = nbformat.read(f, as_version=4) #read the notebook to variable

    # Get a list of Parameter objects
    orig_parameters = extract_parameters(nb)
    
    # Update one or more parameters
    params = parameter_values(orig_parameters, adata_unintegrated= adata_unintegrated, integrated_adata=integrated_adata, tools=tools, data_dir=data_dir, batch=batch, label=label )

    # Make a notebook object with these definitions
    new_nb = replace_definitions(nb, params)

    #Write the edited notebook to file
    nbformat.write(new_nb, fp=os.getcwd() + "/integration_results.ipynb")


#This if statement is executed when script is called
if __name__ == "__main__":

    #Read input arguments
    data_dir = sys.argv[1]
    benchmarking_settings_file = sys.argv[2]
    tools = sys.argv[3:len(sys.argv)] 
    
    
    #As nextflow messes the python list created in parse_methods process, remove the unwanted '[', ',' and ']' -characters from the list items, and update the list
    tools[:] = [re.sub(r'[\W*]', '', tool) for tool in tools]
    integrated_adata = {"scanorama":f"{data_dir}/scanorama_adata.h5ad", "scvi":f"{data_dir}/scvi_adata.h5ad", "combat":f"{data_dir}/combat_adata.h5ad", 
                        "scanvi":f"{data_dir}/scanvi_adata.h5ad", "harmony":f"{data_dir}/harmony_adata.h5ad", "desc":f"{data_dir}/desc_adata.h5ad", 
                        "liger":f"{data_dir}/liger_adata.h5ad", "bbknn":f"{data_dir}/bbknn_adata.h5ad", "scgen":f"{data_dir}/scgen_adata.h5ad",
                        "trvae":f"{data_dir}/trvae_adata.h5ad", "mnn":f"{data_dir}/mnn_adata.h5ad", "saucie":f"{data_dir}/saucie_adata.h5ad",
                        "seurat_cca":f"{data_dir}/seurat_cca_output.rds", "seurat_rpca":f"{data_dir}/seurat_rpca_output.rds",
                        "fastmnn":f"{data_dir}/fastmnn_output.rds", "conos":f"{data_dir}/conos_adata.h5"}
    print(tools)
    print(integrated_adata)
    
    #Read the benchmarking_settings.csv file to pandas df
    settings_df = pd.read_csv(benchmarking_settings_file)
    for tool in tools: 
        method, batch, num_hvg, r_output, seurat_output, scale, label_key = extract_preprocessing_settings(tool, settings_df)  #Extract settings from the input df
        del method, num_hvg, r_output, seurat_output, scale
        break
    
    update_notebook_parameters(f"{data_dir}/unintegrated.h5ad", integrated_adata, tools, data_dir, batch, label_key)
    
   