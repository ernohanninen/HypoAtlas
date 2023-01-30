import scanpy as sc
import numpy as np


def scanorama(data_path, output_path):
    import scanorama
    adata = sc.read(data_path)
    # List of adata per batch
    
    
    batch_cats = adata.obs.batch.cat.categories
    adata_list = [adata[adata.obs.batch == b].copy() for b in batch_cats]
    scanorama.integrate_scanpy(adata_list)

    adata.obsm["scanorama"] = np.zeros((adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))
    for i, b in enumerate(batch_cats):
        adata.obsm["scanorama"][adata.obs.batch == b] = adata_list[i].obsm["X_scanorama"]    
    
    adata.write(output_path + "scanorama_adata.h5ad")

def liger(data_path, output_path):
    
    adata = sc.read(data_path)
    batch_cats = adata.obs.batch.cat.categories
    
    import pyliger
    
    bdata = adata.copy()
    # Pyliger normalizes by library size with a size factor of 1
    # So here we give it the count data
    bdata.X = bdata.layers["counts"]
    # List of adata per batch
    adata_list = [bdata[bdata.obs.batch == b].copy() for b in batch_cats]
    for i, ad in enumerate(adata_list):
        ad.uns["sample_name"] = batch_cats[i]
        # Hack to make sure each method uses the same genes
        ad.uns['var_gene_idx'] = np.arange(bdata.n_vars)


    liger_data = pyliger.create_liger(adata_list, remove_missing=False, make_sparse=False)
    # Hack to make sure each method uses the same genes
    liger_data.var_genes = bdata.var_names
    pyliger.normalize(liger_data)
    pyliger.scale_not_center(liger_data)
    pyliger.optimize_ALS(liger_data, k=30)
    pyliger.quantile_norm(liger_data)


    adata.obsm["liger"] = np.zeros((adata.shape[0], liger_data.adata_list[0].obsm["H_norm"].shape[1]))
    for i, b in enumerate(batch_cats):
        adata.obsm["liger"][adata.obs.batch == b] = liger_data.adata_list[i].obsm["H_norm"]
        
    adata.write(output_path + "liger_adata.h5ad")
    

def harmony(data_path, output_path, batch):
    from harmony import harmonize
    
    adata = sc.read(data_path)
    sc.tl.pca(adata)
    adata.obsm["harmony"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key = batch)
    
    adata.write(output_path + "harmony_adata.h5ad")
    
    

def scvi(data_path, output_path, batch):
    
    import scvi
    
    adata = sc.read(data_path)

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch)
    vae = scvi.model.SCVI(adata, gene_likelihood="nb", n_layers=2, n_latent=30)
    vae.train(max_epochs=10)
    adata.obsm["scvi"] = vae.get_latent_representation()    
    adata.write(output_path + "scvi_adata.h5ad")
    
    
def scanvi(data_path, output_path, batch, label_key):
    import scvi
    
    adata = sc.read(data_path)
    
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch)
    vae = scvi.model.SCVI(adata, gene_likelihood="nb", n_layers=2, n_latent=30)
    vae.train(max_epochs=10)
    
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key=label_key,
        unlabeled_category="Unknown",
    )
    lvae.train(max_epochs=20, n_samples_per_label=100)
    adata.obsm["scanvi"] = lvae.get_latent_representation()
    
    adata.write(output_path + "scanvi_adata.h5ad")
    
    
def desc(data_path, output_path, batch):
    
    #ImportError: You must install pydot (`pip install pydot`) and install graphviz (see instructions at https://graphviz.gitlab.io/download/) for plot_model to work.
    
    import desc, tempfile

    adata = sc.read(data_path)


    temp_dir = tempfile.TemporaryDirectory()
    tmp_dir = temp_dir.name

    adata_out = adata.copy()

    adata_out = desc.scale_bygroup(adata_out, groupby=batch, max_value=6)

    adata_out = desc.train(
        adata_out,
        dims=[adata.shape[1], 128, 32],
        tol=0.001,
        n_neighbors=10,
        batch_size=256,
        louvain_resolution=0.8,
        save_encoder_weights=False,
        save_dir=tmp_dir,
        do_tsne=False,
        use_GPU=False,
        num_Cores=5,
        use_ae_weights=False,
        do_umap=False,
    )

    adata.obsm["desc"] = adata_out.obsm["X_Embeded_z" + str(0.8)]
    
    adata.write(output_path + "desc_adata.h5ad")


def combat(data_path, output_path, batch):
    adata = sc.read(data_path)
    adata_int = adata.copy()
    sc.pp.combat(adata_int, key=batch) # Replaces the adata.X with corrected data
    print(adata_int)
    adata.obsm["combat"] = adata_int.X
    
    adata.write(output_path + "combat_adata.h5ad")


def bbknn(data_path, output_path, batch):
    import bbknn
    
    
    adata = sc.read(data_path)
    sc.tl.pca(adata)
    
    if adata.n_obs < 1e5:
        adata = bbknn.bbknn(adata, batch_key=batch, copy=True)
    if adata.n_obs >= 1e5:
        adata = bbknn.bbknn(adata, batch_key=batch, neighbors_within_batch=25, copy=True)
    
    print("bbknn adata")
    print(adata)
    adata.write(output_path + "bbknn_adata.h5ad")
    

def scgen(data_path, output_path, batch, label_key):
    import scgen 
    
    adata = sc.read(data_path)
    
    scgen.SCGEN.setup_anndata(adata, batch_key=batch, labels_key=label_key)
    
    #Creating the model
    model = scgen.SCGEN(adata)
    
    #Training
    model.train(
        max_epochs=10,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    
    corrected_adata = model.batch_removal()
    print(corrected_adata)
    adata.obsm["scgen"] = corrected_adata.X
    print(adata.obsm["scgen"])
    
    adata.write(output_path + "scgen_adata.h5ad")

def trvae(data_path, output_path, batch, label):
    from scarches.dataset.trvae.data_handling import remove_sparsity
    import torch
    import scarches as sca
    
    adata = sc.read(data_path)
    
    #Get target condition (the dominant batch covariate)
    target_condition = list(adata.obs[batch].value_counts().idxmax())
    
    trvae_epochs = 50
    surgery_epochs = 50

    early_stopping_kwargs = {
        "early_stopping_metric": "val_unweighted_loss",
        "threshold": 0,
        "patience": 20,
        "reduce_lr": True,
        "lr_patience": 13,
        "lr_factor": 0.1,
    }
    
    #Split the dataset int reference and query dataset

    #Process the data for trVAE
    #adata = adata.raw.to_adata()
    adata = remove_sparsity(adata) #if adata.X is sparse matrix -> converts it in to normal matrix

    source_adata = adata[~adata.obs[batch].isin(target_condition)]
    target_adata = adata[adata.obs[batch].isin(target_condition)]

    #Get source conditions (all batches of the data)
    source_conditions = source_adata.obs[batch].unique().tolist()
    
    print("Source adata: ", source_adata)
    print("Target adata: ", target_adata)
    print("Source conditions: ",source_conditions)
    print("Target conditions: ",target_condition)
    
    #Create the TRVAE model
    trvae = sca.models.TRVAE(
        source_adata,
        condition_key=batch,
        conditions=source_conditions,
        hidden_layer_sizes=[128, 128],
    )
    
    #Training trVAE with the reference dataset (source_adata)
    trvae.train(
        n_epochs=trvae_epochs,  
        alpha_epoch_anneal=200,
        early_stopping_kwargs=early_stopping_kwargs,
        seed = 42
    )
    
    adata_latent = sc.AnnData(trvae.get_latent())
    adata_latent.obs[label] = source_adata.obs[label].tolist()
    adata_latent.obs[batch] = source_adata.obs[batch].tolist()
         
       
    print(adata_latent)
    
    print("Fine tune te reference model with query data")
    #Fine tune te reference model with query data
    new_trvae = sca.models.TRVAE.load_query_data(adata=target_adata, reference_model=trvae)
    
    new_trvae.train(n_epochs=surgery_epochs,alpha_epoch_anneal=200,early_stopping_kwargs=early_stopping_kwargs,weight_decay=0, seed = 42)
    
    adata_latent = sc.AnnData(new_trvae.get_latent())
    adata_latent.obs[label] = target_adata.obs[label].tolist()
    adata_latent.obs[batch] = target_adata.obs[batch].tolist()
    
    
    # Get latent representation of reference + query dataset 
    full_latent = sc.AnnData(new_trvae.get_latent(adata.X, adata.obs[batch]))
    full_latent.obs[label] = adata.obs[label].tolist()
    full_latent.obs[batch] = adata.obs[batch].tolist()
    
    print(full_latent)
    
    adata.obsm["trvae"] = full_latent.X
    
    adata.write(output_path + "trvae_adata.h5ad")
    
        
def mnn(data_path, output_path, batch, hvg=None):
    import mnnpy
    adata = sc.read(data_path)
    print("Before integration")
    print(adata)
    
    #Split batches for mnn algorithm
    split = []
    batch_categories = adata.obs[batch].cat.categories
    if hvg is not None:
        adata = adata[:, hvg]
    for i in batch_categories:
        split.append(adata[adata.obs[batch] == i].copy())
    
    print(batch_categories)
    print(split)
    
    corrected, _, _ = mnnpy.mnn_correct(
        *split,
        var_subset=hvg,
        batch_key=batch,
        batch_categories=batch_categories,
        index_unique=None
    )
    print("After integration")
    print(corrected)
    #Print the countes in adata.X
    adata.X.todense()[185:190,185:190]
    corrected.X.todense()[185:190,185:190]
    
    adata.obsm["mnn"] = corrected.X
    adata.write(output_path + "mnn_adata.h5ad")
    
    
    
def saucie(data_path, output_path, batch):

    adata = sc.read(data_path)
    import sys
    sys.path.insert(0, "../../../../") #Path in where saucie algorithm is cloned, saucie is imported from this path
    import SAUCIE
    import sklearn.decomposition
    import scipy as sp
    from scipy.sparse import issparse

    pca_op = sklearn.decomposition.PCA(100)
    if isinstance(adata.X, sp.sparse.csr_matrix):
        expr = adata.X.A
    else:
        expr = adata.X
        
    data = pca_op.fit_transform(expr)
    saucie = SAUCIE.SAUCIE(100, lambda_b=0.1)
    loader_train = SAUCIE.Loader(data, labels=adata.obs[batch].cat.codes, shuffle=True)
    loader_eval = SAUCIE.Loader(data, labels=adata.obs[batch].cat.codes, shuffle=False)
    saucie.train(loader_train, steps=5000)
    ret = adata.copy()
    ret.obsm["X_emb"] = saucie.get_reconstruction(loader_eval)[0]
    adata.obsm["saucie"] = pca_op.inverse_transform(ret.obsm["X_emb"])

    adata.write(output_path + "saucie_adata.h5ad")


    
    
    
    