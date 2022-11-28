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

#Function which extracts the preprocessing settings for each tool
def extract_preprocessing_settings(tool, settings_df):
    #Get the settings of specific tool and read them to variable
    method, batch, num_hvg, r_output, seurat_output, scale, label_key, epochs = settings_df.loc[settings_df["integration_method"] == tool].values[0]
    
    num_hvg = int(num_hvg)

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
    print(method, batch, num_hvg, r_output, seurat_output, scale)
    return  method, batch, num_hvg, r_output, seurat_output, scale #return settings

#Function which updates notebook parameters
def update_notebook_parameters(dataset, integration_method, batch, num_hvg, r_output, seurat_output, scale):
    #Open the preprocessing notebook
    with open("../../../BenchmarkStudy/Notebooks_to_run/run_preprocessing.ipynb") as f:
        nb = nbformat.read(f, as_version=4) #read the notebook to variable

    # Get a list of Parameter objects
    orig_parameters = extract_parameters(nb)
    

    # Update one or more parameters
    params = parameter_values(orig_parameters, method=integration_method,  file="/home/ernohanninen/master_project/Data/lungatlas.h5ad",  out = os.getcwd(), hvg=num_hvg, rout=r_output, batch=batch, seurat = seurat_output, scale = scale)

    # Make a notebook object with these definitions
    new_nb = replace_definitions(nb, params)

    #Write the edited notebook to file
    nbformat.write(new_nb, fp="../../../BenchmarkStudy/Notebooks_to_run/run_preprocessing.ipynb")
    

#Function which performs data preprocessing
def run_preprocessing(integration_method, r_output):
    #using jupyter nbconvert bash-command to execute the notebook, store the executed notebook to working dir 
    bashCommand = "jupyter nbconvert --to notebook --execute ../../../BenchmarkStudy/Notebooks_to_run/run_preprocessing.ipynb --output " + os.getcwd() + "/" + integration_method + "_run_preprocessing.ipynb"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE) #Run command
    output, error = process.communicate() #Get the output
    print("OUTPUT : ", output)
    print(type(r_output))
    print(r_output)
    if r_output == True:
        return os.getcwd() + "/" + integration_method + "_seurat.rds"
    return os.getcwd() + "/" + integration_method + "_adata.h5ad"
    

#This if statement is executed when script is called
if __name__ == "__main__":

    #Read input arguments
    dataset = sys.argv[1]
    benchmarking_settings_file = sys.argv[2]
    tools_to_benchmark = sys.argv[3:int(len(sys.argv))]
    
    #Read the benchmarking_settings.csv file to pandas df
    settings_df = pd.read_csv(benchmarking_settings_file)
    
    adata_dict = {} #store the processed adata file paths to dict
    
    #loop over the tools, during every iteration call the extract_preprocessing_settings and run_preprocessing functions
    for tool in tools_to_benchmark: 
        tool = re.sub(r'[\W_]', '', tool) #As nextflow messes the python list created in parse_methods process, remove the unwanted '[', ',' and ']' -characters from the list items    
        method, batch, num_hvg, r_output, seurat_output, scale = extract_preprocessing_settings(tool, settings_df)  
        update_notebook_parameters(dataset, method, batch, num_hvg, r_output, seurat_output, scale)
        file_path = run_preprocessing(method, r_output)
        adata_dict.update({tool:file_path}) #Store the preprocessed adata file path to dict 
    
    #Store the file paths to csv file
    with open("preprocessed_adata.csv", "w") as csvfile:
        fieldnames= ["integration_method", "adata_path"] #Initialize the headers
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames) #Create writer object
        writer.writeheader()
        for key, value in adata_dict.items(): #Loop over the dict and store the key and value to file
            writer.writerow({"integration_method":key, "adata_path":value})    