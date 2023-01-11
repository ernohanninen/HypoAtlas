#!/usr/bin/env python
"""
Author: Erno Hänninen
Created: 2022-12-02
title collect_integration_metrics.py

Description:
- 
"""

import csv
import sys, re
import pandas as pd

data, num_columns = [], -1 #initialize variables. Data is used as list of lists in where each list contains metrics for from integration algorithm. num_columns is used as global variable

#This function is callled once per each benchmarked integration algorithm
def extract_data(filename): 
    
    #Read the integration metrics to dataframe
    integration_df = pd.read_csv("../../../BenchmarkStudy/Output/Integration/Metrics/"+filename)
    
    global num_columns #num_columns global variable is used to keep in track how many metric columns there where in the csv file for the previous integration method
    #If the amount of columns doesn't match, print error message and stop the program. This is done to avoid error when merging all the data from different files to one csv file
    if num_columns != len(list(integration_df.columns)) and num_columns != -1:
        print("The number of columns in csv files doesn't match, the program is terminated...")
        sys.exit()
    num_columns = len(list(integration_df.columns)) #Update the num_columns global variable
    
    df_columns = list(integration_df.columns) #Convert the dataframe columns to list
    df_columns.insert(0, "Method") #Add method column to the list
    
    return list(integration_df.values[0]), df_columns #Return the values of df converted to list and the df_columns list
    

#Function that appends the data list of lists, called during every iteration
def collect_data(tool, integration_metrics, df_columns):
    integration_metrics.insert(0, tool) #Add the integration method the integration_metrics list
    data.append(integration_metrics) #Update the data list with list
    return df_columns


#Function that creates pandas dataframe from values stored in data list of lists and columns stored in df_columns list 
def create_df(data, df_columns):
    df = pd.DataFrame(data, columns = df_columns) #Create the df
    
    #Write the df to csv file. Ignore the row index and fill nan values with NA
    df.to_csv("merged_benchmarking_metrics.csv", index=False, na_rep="NA")


#This if statement is executed when script is called
if __name__ == "__main__":
    
    tools = sys.argv[1:len(sys.argv)] #read the input arguments
    
    #As nextflow messes the python list created in parse_methods process, remove the unwanted '[', ',' and ']' -characters from the list items, and update the list
    tools[:] = [re.sub(r'[\W*]', '', tool) for tool in tools]
    
    for tool in tools: # loop over the tools list and call the functions during every iteration
        extracted_data =  extract_data(tool + "_integration_metrics.csv")
        df_columns = collect_data(tool, extracted_data[0], extracted_data[1])
    
    create_df(data, df_columns) #Function that uses the data from extract_data and collect_data functions to construnct the df
    