#!/usr/bin/env python
"""
Author: Erno Hänninen
Created: 2022-11-28
title jupyter_functions.py

Description:
- 
"""

import os

#cwd = os.getcwd()
#os.chdir("/home/bns631/HypoAtlas")
import scanpy as sc
#import scib
import pandas as pd
import h5py
import scipy.io
from scipy.sparse import csr_matrix

#os.chdir(cwd)

def read_conos_to_scanpy(conos_path, original_adata):
        
    #As the output of saveConosForScanPy wasn't directly readable to scanpy, reading the matrixes of conos_adata.h5 file independently and constructing the adata object from matrixes and dataframes
    dictionary = {}
    with h5py.File(conos_path, "r") as f: #Open file  
        #Read the metadata and extract the cell id's
        metadata = pd.DataFrame(f["metadata"]['metadata.df'][:])
        metadata.index = metadata.CellId
        del metadata["CellId"]
        
        #read genes, embedding and pseudopca matrices as pandas dataframe
        gene_df = pd.DataFrame(f["genes"]["genes.df"][:]) # Creates a df of the returned numpy array
        embedding_df = pd.DataFrame(f["embedding"]["embedding.df"][:])  
        pseudopca_df = pd.DataFrame(f["pseudopca"]["pseudopca.df"][:])

        #Construct the graph connectivity matrix
        shape = (f["graph_connectivities"]["shape"][:][0],f["graph_connectivities"]["shape"][:][1])
        graph_conn_mtx = csr_matrix((f["graph_connectivities"]["data"][:], f["graph_connectivities"]["indices"][:], f["graph_connectivities"]["indptr"][:]), shape=shape)
        
        #Construct the graph distance matrix
        shape = (f["graph_distances"]["shape"][:][0],f["graph_distances"]["shape"][:][1])
        graph_dist_mtx = csr_matrix((f["graph_distances"]["data"][:], f["graph_distances"]["indices"][:], f["graph_distances"]["indptr"][:]), shape=shape)
        
        #Create the count matrix
        shape = (f["raw_count_matrix"]["shape"][:][0],f["raw_count_matrix"]["shape"][:][1])
        count_csr_matrix = csr_matrix((f["raw_count_matrix"]["data"][:], f["raw_count_matrix"]["indices"][:], f["raw_count_matrix"]["indptr"][:])).transpose()
        
        #Store the count matrix to file and then read it and construct an adata object  
        scipy.io.mmwrite("raw_count_matrix.mtx", count_csr_matrix)
        adata = sc.read_mtx("raw_count_matrix.mtx")
        os.remove("raw_count_matrix.mtx") #Remove the created temp file
        
        #Initialize var and obs columns
        adata.var_names = gene_df["gene"].values
        adata.obs_names = metadata.index.values
        adata.obs = metadata.copy()
        #Pca column
        adata.X_pca = pseudopca_df.values
        adata.obsm["X_pca"] = pseudopca_df.values
        #Umap column
        adata.X_umap = embedding_df.values
        adata.obsm["X_umap"] = embedding_df.values
        #Neighbors column
        adata.uns["neighbors"] = dict(
            connectivities=graph_conn_mtx.tocsr(), distances=graph_dist_mtx.tocsr()
        )
        
        #In the converison some of the columns are stored as byte objects, decoding the to UTF-8 format
        adata.var.index = adata.var.index.str.decode("utf-8")
        adata.obs.index = adata.obs.index.str.decode("utf-8")

        #In the conversion some the categorical columns gets numeric value, therfore taking the categorical columns from the original adata to the integrated adata
        for column in original_adata.obs.columns:
            if original_adata.obs[column].dtypes == "category":
                adata.obs[column] = original_adata.obs[column]
                
    return adata

# Function which plots the result of data interation
def plot_results(adata, method, result):

    #Prepare the dataset
    if result == 'embed':
        #adata.obsm["X_emb"] = adata.obsm["X_pca"].copy() # This line is not executed in the scib-pipeline, if the pca is not computed this will fail   
        scib.pp.reduce_data(adata, n_top_genes=None, neighbors=True, use_rep='X_emb', pca=False, umap=False)
    elif result == 'full':
        sc.pp.filter_genes(adata, min_cells=1)
        scib.pp.reduce_data(adata, n_top_genes=2000, neighbors=True, use_rep='X_pca', pca=True, umap=False)
        
    # Calculate embedding
    if method == 'conos':
        print('Calculating graph embedding...')
        sc.tl.draw_graph(adata, key_added_ext='graph')
        basis = 'draw_graph_graph'
        label = 'Graph'
    else: #Requires neighbors to be computed
        print('Calculating UMAP...')
        sc.tl.umap(adata)
        basis = 'umap'
        label = 'UMAP'
        
    return adata, basis
        
    
#Function which runs the data integration metrics and stores the result to file
def compute_metrics(method, adata, adata_int, batch_key, label_key, embed, integration_type):
    if label_key in adata.var_names or label_key in adata.obs.columns:
        #Doesn't contain cell cycle metrictraject
        cluster_key = label_key
        metrics_df = scib.metrics.metrics(adata, adata_int, batch_key, label_key, embed='X_pca', cluster_key=cluster_key, cluster_nmi=None, ari_=True, nmi_=True, nmi_method='arithmetic', nmi_dir=None, silhouette_=True, si_metric='euclidean', pcr_=True, cell_cycle_=False, organism='human', hvg_score_=False, isolated_labels_=True, isolated_labels_f1_=True, isolated_labels_asw_=True, n_isolated=None, graph_conn_=True, trajectory_=False, kBET_=False, lisi_graph_=True, ilisi_=True, clisi_=True, subsample=0.5, n_cores=1, type_="full", verbose=False)
        #metrics_df = scib.metrics.metrics(adata, adata_int, batch_key, label_key, embed='X_pca', cluster_key=cluster_key, cluster_nmi=None, ari_=False, nmi_=False, nmi_method='arithmetic', nmi_dir=None, silhouette_=False, si_metric='euclidean', pcr_=False, cell_cycle_=False, organism='human', hvg_score_=False, isolated_labels_=False, isolated_labels_f1_=False, isolated_labels_asw_=False, n_isolated=None, graph_conn_=False, trajectory_=False, kBET_=False, lisi_graph_=False, ilisi_=False, clisi_=False, subsample=0.5, n_cores=1, type_="full", verbose=False)
        
        metrics_df.dropna(inplace=True)
        metrics_df = metrics_df.transpose()
        print(metrics_df)  
    
    #If dataset doesn't contain cell type labels compute metrics that doesn't require label information
    else:
        hvg_overlap_value = scib.metrics.hvg_overlap(adata, adata_int, batch_key, n_hvg=200)
        pcr_score = scib.metrics.pcr_comparison(adata,adata_int,covariate=batch_key,embed=embed,n_comps=50)
        ilisi_score = scib.metrics.ilisi_graph(adata_int, batch_key, k0=90, type_=integration_type, subsample=None, scale=True, n_cores=1, verbose=False)
        
        #Round the values, store them to dictinary, and use the dictionary to construct dataframe
        metrics_dict = {"hvg_overlap_value":[round(hvg_overlap_value,4),], "pcr_score":[round(pcr_score,4),],"ilisi_score":[round(ilisi_score,4),]}        
        metrics_df = pd.DataFrame(metrics_dict)
    
    #Store dataframe to csv file
    metrics_df.to_csv("../Metrics/"+ method + "_integration_metrics.csv", index=False)
    return metrics_df #return df

