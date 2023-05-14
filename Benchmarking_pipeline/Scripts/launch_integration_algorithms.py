"""
Author: Erno Hänninen

Created: 14.01.2023

Title: launch_integration_algorithms.py

Description:
- Calls the integration algorithms written in python from integration_algorithms.py file. 

Procedure
- Takes arguments from nextflow pipeline
- Reads the integration algorithm parameters from the configuration file
- Executes the desired integration algorithm by calling it from integration_algorithms.py script
    
List of non-standard modules:
- pandas

Usage:
- This script is launched from the pipeline (processes.nf)
"""

from integration_algorithms import *
import sys, csv
from run_preprocessing import extract_preprocessing_settings
import pandas as pd

# Function that calls algorithms from integration_algorithms.py
def call_algorithm(data_path, tool, output_path, batch, label_key):
    if tool == "scanorama": scanorama(data_path, output_path,batch)
    elif tool == "liger": liger(data_path, output_path,batch)
    elif tool == "harmony": harmony(data_path, output_path, batch)
    elif tool == "scvi": scvi(data_path, output_path, batch)
    elif tool == "scanvi": scanvi(data_path, output_path, batch, label_key)
    elif tool == "desc": desc(data_path, output_path, batch)
    elif tool == "combat": combat(data_path, output_path, batch)
    elif tool == "bbknn": bbknn(data_path, output_path, batch)
    elif tool == "scgen": scgen(data_path, output_path, batch, label_key)
    elif tool == "trvae": trvae(data_path, output_path, batch, label_key)
    elif tool == "mnn": mnn(data_path, output_path, batch)
    elif tool == "saucie": saucie(data_path, output_path, batch)
        
# This if statement is executed when script is launched
if __name__ == "__main__":
    
    # Read input arguments
    preprocessed_adata_files = sys.argv[1]
    tool = sys.argv[2]
    output_path = sys.argv[3]
    benchmarking_settings_file = sys.argv[4]
     
    # Read the paths preprocessed data from file
    with open(preprocessed_adata_files, "r") as csvfile:
        reader = csv.reader(csvfile) #Initialize csv reader
        next(reader) #Read the header
        for row in reader: #Loop thru the rows in the file
            if row[0] == tool:
                data_path = row[1]
                break
    
    # Read the benchmarking_settings.csv file to pandas df
    settings_df = pd.read_csv(benchmarking_settings_file)

    # Extract the batch and cell type information from the configuration file
    method, batch, num_hvg, r_output, seurat_output, scale, label_key = extract_preprocessing_settings(tool, settings_df)  #Extract settings from the input df
    del method, num_hvg, r_output, seurat_output, scale

    call_algorithm(data_path, tool, output_path, batch, label_key)
