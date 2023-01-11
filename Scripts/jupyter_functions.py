#!/usr/bin/env python
"""
Author: Erno Hänninen
Created: 2022-11-28
title jupyter_functions.py

Description:
- 
"""

import os

cwd = os.getcwd()
os.chdir("/home/ernohanninen/master_project")
import scanpy as sc
import scib
import pandas as pd
os.chdir(cwd)

# Function which plots the result of data interation
def plot_results(adata, method, result):

    #Prepare the dataset
    if result == 'embed':
        adata.obsm["X_emb"] = adata.obsm["X_pca"].copy() # This line is not executed in the scib-pipeline, if the pca is not computed this will fail   
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

