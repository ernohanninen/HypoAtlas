

import scanpy as sc
import scib
import subprocess

def activate_env(env):
    subprocess.run('conda activate ' + env, shell=True)
    import os
    print(os.environ['CONDA_DEFAULT_ENV'])

def deactivate_env():
    subprocess.run('source deactivate', shell=True)


def plot_results(adata, method, result):
    
    adata = sc.read_h5ad(adata)
    
    #Prepare the dataset
    if result == 'embed':
        scib.pp.reduce_data(adata, n_top_genes=None, neighbors=True, use_rep='X_emb', pca=False, umap=False)
    elif result == 'full':
        sc.pp.filter_genes(adata, min_cells=1)
        scib.pp.reduce_data(adata, n_top_genes=2000, neighbors=True, use_rep='X_pca', pca=True, umap=False)
        
    
    # Calculate embedding
    if method.startswith('conos'):
        print('Calculating graph embedding...')
        sc.tl.draw_graph(adata, key_added_ext='graph')
        basis = 'draw_graph_graph'
        label = 'Graph'
    else:
        print('Calculating UMAP...')
        sc.tl.umap(adata)
        basis = 'umap'
        label = 'UMAP'
    
    return adata, basis
        
    
    
    
def compute_metrics():
    print("OK")