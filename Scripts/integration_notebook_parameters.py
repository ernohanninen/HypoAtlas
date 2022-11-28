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



#Function which extracts the settings for the integration tool
def extract_benchmarking_settings(tool, settings_df):
    #Get the settings of specific tool and read them to variable
    method, batch, num_hvg, r_output, seurat_output, scale,label_key,epochs = settings_df.loc[settings_df["integration_method"] == tool].values[0]
    print("SETTINGS TO BE USED:")
    
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
    #Open the preprocessing notebook
    with open("../../../BenchmarkStudy/Notebooks_to_run/run_" + integration_method + ".ipynb") as f:
        nb = nbformat.read(f, as_version=4) #read the notebook to variable

    # Get a list of Parameter objects
    orig_parameters = extract_parameters(nb)
    
    print("SETTINGS: ", integration_method, adata_path, batch_key, label_key, num_hvg, epochs)
    # Update parameters
    if integration_method == "scvi":
        params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, labels=label_key, epochs=epochs)
    elif integration_method == "scanvi":
        params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, labels=label_key, epochs=epochs)
    elif integration_method == "scgen":
        params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, labels=label_key, epochs=epochs)
    elif integration_method == "scanorama":
        params = parameter_values(orig_parameters, input_data_path=adata_path,  batch=batch_key, labels=label_key)
        
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
    
    with open(preprocessed_adata_files, "r") as csvfile:
        reader = csv.reader(csvfile) #Initialize csv reader
        next(reader) #Read the header
        for row in reader: #Loop thru the rows in the file
            print(row)
            method,batch,label_key,num_hvg,epochs = extract_benchmarking_settings(row[0], settings_df)
            update_notebook_parameters(method, row[1], batch, label_key, num_hvg, epochs)    
   
   
   