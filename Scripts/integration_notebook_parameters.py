#!/usr/bin/env python
"""
Author: Erno Hänninen
Created: 2022-11-25
title run_integration_notebooks.py

Description:
- 
"""

import os
print(os.environ['CONDA_DEFAULT_ENV'])

import sys, csv
import subprocess #Allows running bash
import os
import pandas as pd
import re
#packages which allows edit notebooks from python
import nbclient
import nbformat
from nbparameterise import (
    extract_parameters, replace_definitions, parameter_values
)
import papermill as pm



#Function which extracts the settings for the integration tool
def extract_benchmarking_settings(tool, settings_df):
    #Get the settings of specific tool and read them to variable
    method, batch, num_hvg, r_output, seurat_output, scale,label_key,epochs = settings_df.loc[settings_df["integration_method"] == tool].values[0]
    print("SETTINGS TO BE USED:")
    
    #If statements that converts the string true and false values to boolean, and numeric string values to integer
    if epochs == "None":
        epochs = None
    else:
        epochs = int(epochs)
    if r_output == 'True':
        r_output = True
    elif r_output == 'False':
        r_output = False
    
    print(method,batch,label_key,num_hvg,epochs)
    return  method, batch,label_key,num_hvg,epochs #return settings

#Function which updates notebook parameters
def update_notebook_parameters(integration_method, adata_path, batch_key, label_key, num_hvg, epochs):
    
    #For notebooks written with R using papermill python package to parameterize the notebooks. The papermill requires that the variables which takes values from papermill are tagged in the notebook
    #Also the downstream quantification notebooks of these integration algorithms are parameterized using papermill, event these are notebooks written in python
    if integration_method == "conos":
        # Integration notebook
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb", os.getcwd() + "/run_" + integration_method + ".ipynb", parameters = dict(input_data_path=adata_path, hvg=adata_path+"_hvg.RDS", batch=batch_key, label=label_key), prepare_only=True, language="R")
        # quantification notebook  
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/quantify_conos.ipynb", os.getcwd() + "/quantify_conos.ipynb", parameters = dict(batch=batch_key, label=label_key, original_data=adata_path), prepare_only=True)
        
    elif integration_method == "fastmnn":
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb", os.getcwd() + "/run_" + integration_method + ".ipynb", parameters = dict(input_data_path=adata_path, hvg=adata_path+"_hvg.RDS", batch=batch_key, label=label_key), prepare_only=True, language="R")
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/quantify_fastmnn.ipynb", os.getcwd() + "/quantify_fastmnn.ipynb", parameters = dict(batch=batch_key, label=label_key, original_data=adata_path), prepare_only=True)
        
    elif integration_method == "liger":
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb", os.getcwd() + "/run_" + integration_method + ".ipynb", parameters = dict(input_data_path=adata_path, hvg=adata_path+"_hvg.RDS", batch=batch_key, label=label_key), prepare_only=True, language="R")
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/quantify_liger.ipynb", os.getcwd() + "/quantify_liger.ipynb", parameters = dict(batch=batch_key, label=label_key, original_data=adata_path), prepare_only=True)
        
    elif integration_method == "seurat_cca":
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb", os.getcwd() + "/run_" + integration_method + ".ipynb", parameters = dict(input_data_path=adata_path, hvg=adata_path+"_hvg.RDS", batch=batch_key, label=label_key), prepare_only=True, language="R")
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/quantify_seurat_cca.ipynb", os.getcwd() + "/quantify_seurat_cca.ipynb", parameters = dict(batch=batch_key, label=label_key, original_data=adata_path), prepare_only=True)
        
    elif integration_method == "seurat_rpca":
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb", os.getcwd() + "/run_" + integration_method + ".ipynb", parameters = dict(input_data_path=adata_path, hvg=adata_path+"_hvg.RDS", batch=batch_key, label=label_key), prepare_only=True, language="R")
        pm.execute_notebook("../../../BenchmarkStudy/Notebooks_to_run/quantify_seurat_rpca.ipynb", os.getcwd() + "/quantify_seurat_rpca.ipynb", parameters = dict(batch=batch_key, label=label_key, original_data=adata_path), prepare_only=True)
    
    # THe integration notebooks written in python are parameterized using nbparameterise package
    else:
        #Open the preprocessing notebook
        with open("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb") as f:
            nb = nbformat.read(f, as_version=4) #read the notebook to variable

        # Get a list of Parameter objects
        orig_parameters = extract_parameters(nb)
        
        print("SETTINGS: ", integration_method, adata_path, batch_key, label_key, num_hvg, epochs)
        # Update parameters, different alogirhms needs different parameters, if statements takes care which parameters to use
        if integration_method == "scvi" or integration_method == "scanvi" or integration_method == "scgen":
            params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, label=label_key, epochs=epochs)
            
        elif integration_method == "desc":
            params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, label=label_key, use_gpu = False, ncores = 1)   
            
        elif integration_method == "scanorama" or integration_method == "combat" or integration_method == "harmony" or integration_method == "mnn" or integration_method == "saucie" or integration_method == "trvae" or integration_method == "trvae_tl" or integration_method == "bbknn":
            params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, label=label_key)

        # Make a notebook object with these definitions
        new_nb = replace_definitions(nb, params)

        #Write the edited notebook to file
        nbformat.write(new_nb, fp=os.getcwd() + "/run_" + integration_method + ".ipynb")
    

#This if statement is executed when script is called
if __name__ == "__main__":

    #Read input arguments
    preprocessed_adata_files = sys.argv[1]
    benchmarking_settings_file = sys.argv[2]
    
    #Read the benchmarking_settings.csv file to pandas df
    settings_df = pd.read_csv(benchmarking_settings_file)
    
    #preprocessed_adata_files csv file contains paths to adata / seurat objects preprocessed for different integration methods the user wants to benchmark
    with open(preprocessed_adata_files, "r") as csvfile:
        reader = csv.reader(csvfile) #Initialize csv reader
        next(reader) #Read the header
        for row in reader: #Loop thru the rows in the file
            method,batch,label_key,num_hvg,epochs = extract_benchmarking_settings(row[0], settings_df) #Extract the data from settings_df dataframe, row[0] argument contains the integration tool
            update_notebook_parameters(method, row[1], batch, label_key, num_hvg, epochs) # Function which updates the notebook parameters for each integration algorithm to be benchmarked