from integration_algorithms import *
import sys, csv
from run_preprocessing import extract_preprocessing_settings
import pandas as pd



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
    
    
    

if __name__ == "__main__":
    
    #Read input arguments
    preprocessed_adata_files = sys.argv[1]
    tool = sys.argv[2]
    output_path = sys.argv[3]
    benchmarking_settings_file = sys.argv[4]
    
    
    with open(preprocessed_adata_files, "r") as csvfile:
        reader = csv.reader(csvfile) #Initialize csv reader
        next(reader) #Read the header
        for row in reader: #Loop thru the rows in the file
            if row[0] == tool:
                data_path = row[1]
                break
    
    #Read the benchmarking_settings.csv file to pandas df
    settings_df = pd.read_csv(benchmarking_settings_file)

    method, batch, num_hvg, r_output, seurat_output, scale, label_key = extract_preprocessing_settings(tool, settings_df)  #Extract settings from the input df
    del method, num_hvg, r_output, seurat_output, scale


    call_algorithm(data_path, tool, output_path, batch, label_key)
