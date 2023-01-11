# Master project


## Data integration method benchmarking study

Benchmarking different data integration methods.

### Benchmark study workflow

Describing the workflow of the benchmarking study.

#### Installations


##### Python environments


In the benchmarking study multiple Conda environments were used to handle the package dependencies. Due to installation issues and package conflicts, each integration algorithm has its own Conda environment, as suggested in [scib documentation](https://scib.readthedocs.io/en/latest/installation.html#installing-additional-packages). However scvi, scanvi, scgen, sctrvaep, and ComBat integration algorithms makes the expection as their dependencies are located in PYenv. In addition to the benchmarking, PYenv, the main environment of the pipeline, is used to manage dependencies in following Nextflow processes used in the benchmarking study: parse_methods, preprocessing, and integration_parameters. 

The Nextflow pipeline specifies the conda environment used to manage the dependencies of nextflow process.

Below the bash commands used to build the different environments, as well as, the environment used to manage the dependencies of the integration algorithm is specified. 


To build the Conda PYenv environment following command was executed:
`conda create --name PYenv pip && conda activate PYenv && conda install pip pymde==0.1.18 scvi-tools=0.19.0 jupyter==1.0.0 nbconvert==7.2.5 && pip install nbparameterise==0.5 scib==1.0.4 && conda install -c bioconda bioconductor-s4vectors && conda install -c bioconda bioconductor-singlecellexperiment==1.20.0 && conda install -c bioconda r-seurat==4.3.0 && pip install scgen trvaep && pip install papermill && pip install scrublet && pip install scArches`



The conda environments used in benchmarking data integration tools:
- ComBat: PYenv, for ComBat algorithm scib uses scanpy implementation. Therefore using PYenv containing SCANPY package to benchmark ComBat


- mnnpy: mnn_env, created by running:
`conda create --name mnn_env jupyter==1.0.0 mnnpy scanpy pip && conda activate mnn_env && pip install scib && conda deactivate`
conda create --name mnn_env jupyter==1.0.0 mnnpy scanpy pip scipy==1.8.1 && conda activate mnn_env && pip install scib && conda deactivate


- desc: desc_env, created by running:
`conda create --name desc_env pip python && conda activate desc_env && pip install desc jupyter tensorflow && git clone https://github.com/theislab/scib.git && cd scib && pip install -e . && conda install rpy2 && conda deactivate`


- saucie: saucie_env, created by running:
1. installing the saucie dependencies, scib and jupyter, and cloning the repo:
`conda create --name saucie_env pip tensorflow==1.14 python==3.7 pandas==1.3.5 && conda activate saucie_env && pip install fcsparser fcswrite numpy scikit-learn matplotlib jupyter && git clone https://github.com/theislab/scib.git && cd scib && pip install -e . && cd .. && git clone https://github.com/KrishnaswamyLab/SAUCIE && pip install louvain && conda deactivate`

- bbknn: bbknn_env, created by running:
`conda create --name bbknn_env jupyter==1.0.0 pip && conda activate bbknn_env && git clone https://github.com/theislab/scib.git && cd scib && pip install -e . && pip install scib[bbknn] && pip install louvain && conda deactivate`


- harmony: harmony_env in where the scib package is installed from source as the scib 1.0.4 version from pypi doesn't contain harmony. The environment was created by running: 
`conda create --name harmony_env jupyter==1.0.0 pip && conda activate harmony_env && git_clone git@github.com:theislab/scib.git && cd scib && pip install . && pip install harmony-pytorch && conda deactivate`


- scanorama: scanorama_env, created by running:
`conda create --name scanorama_env python-annoy pip && conda activate scanorama_env && pip install scib jupyter && git clone https://github.com/brianhie/scanorama.git && cd scanorama/ && python setup.py install --user && conda install R && conda deactivate`


- SCVI: PYenv, scvi-tools are already installed in PYenv, therefore using it to run scvi-algorithm


- SCANVI: PYenv, scvi-tools are already installed in PYenv, therefore using it to run scanvi-algorithm


- scGen: PYenv, scgen is installed in PYenv, therefore using it to run scGen algorithm


- trvae and trvae_tl: scArches, that contains trvae algorithm is installed in PYenv, therefore using it to run trvae algorithm



##### R environment

To build the conda Renv environment follwoing commands was executed:
Set up R Jupyter Notebooks in VS code 
1. Creating R env:
`conda create --name Renv && conda activate Renv && conda install r-recommended r-irkernel jupyter && conda install r-seurat nbconvert r-conos && conda install bioconductor-rhdf5 && conda install -c bioconda bioconductor-batchelor`

2. Add the R-kernel to jupyter by installing kernelspec:
R -e 'IRkernel::installspec()




##### Preparing the locally created and tested conda environments to server:

###### Create first YAML files from the environments:
1. Create a folder for conda environments and navigate there
`mkdir ~/master_project/CondaEnvironments && cd ~/master_project/CondaEnvironments`

2. desc_env:
`conda activate desc_env && conda env export > desc_env.yml && conda deactivate`

3. PYenv:
`conda activate PYenv && conda env export > PYenv.yml && conda deactivate`

4. saucie_env:
`conda activate saucie_env && conda env export > saucie_env.yml && conda deactivate`

5. bbknn_env:
`conda activate bbknn_env && conda env export > bbknn_env.yml && conda deactivate`

6. harmony_env
`conda activate harmony_env && conda env export > harmony_env.yml && conda deactivate`

7. Renv
`conda activate Renv && conda env export > Renv.yml && conda deactivate`

8. scanorama_env
`conda activate scanorama_env && conda env export > scanorama_env.yml && conda deactivate`

9. combat_env
`conda activate combat_env && conda env export > combat_env.yml && conda deactivate`

When the pipeline is cloned to server, the environments can be created from the yaml files.





#### Data
The human hypothalamus data from Herb et al. study was downloaded from https://assets.nemoarchive.org/dat-8ovb8mx

1. Navigate to the folder where the dataset is loaded

2. Unpack the .tgz file:
`tar zxvf Analysis_Herb_Ament_Human_Development_Hypothalamus_Counts.tgz`


2. bdbag package is required to fetch the data from NeMo archieve. Installing the packages required to base conda environment: 
`pip install bdbag && pip install bdbag[boto,globus]`


3. The bdbag is saved to `~/miniconda3/lib/python3.8/site-packages/bdbag`, adding this path to PATH variable  
`export PATH=$PATH:~/miniconda3/lib/python3.8/site-packages/bdbag`


4. Fetch the data:
`bdbag --resolve-fetch all Analysis_Herb_Ament_Human_Development_Hypothalamus_Counts/Raw_data_Herb_Ament_Human_Development_Hypothalamus_Counts`

5. Navigate to the data folder
`cd Analysis_Herb_Ament_Human_Development_Hypothalamus_Counts/Raw_data_Herb_Ament_Human_Development_Hypothalamus_Counts/data`

6. Unzip the hypothalamus files:
`tar -xzf CS13_prosencephalon.mex.tar.gz && tar -xzf CS14_cortex.mex.tar.gz && tar -xzf CS15_forebrain.mex.tar.gz && tar -xzf CS20_hypothalamus1.mex.tar.gz && tar -xzf CS22_2_hypothalamus.mex.tar.gz && tar -xzf CS22_Hypothalamus.mex.tar.gz && tar -xzf GW18_hypothalamus.mex.tar.gz && tar -xzf GW19_hypothalamus.mex.tar.gz && tar -xzf GW20_hypothalamus.mex.tar.gz && tar -xzf GW22T_hypo1.mex.tar.gz && tar -xzf GW25_3V_hypo.mex.tar.gz`



#### Steps

To validate the training process of scgen, scvi, scanvi and trvaep algorithms and to ensure the model achieved convergence, the source code of scib (v. 1.0.4) package was edited:

1. Navigate to the edited file:
`nano ~/miniconda3/envs/PYenv/lib/python3.10/site-packages/scib/integration.py `

2. The original script returns only either integrated adata or the trained model. However we were interested of both as the integrated adata object is needed for plotting and computing the evaluation metrics, and the model is needed for validaton of the training process. Therefore the following changes where done:

 - Return value of scgen() function (row 194): `return corrected_adata, model`
 - Return value of scvi() function (row 246): `return adata, vae`
 - Return value of scanvi() function (row 291) `return adata, scanvae` 


3. In addition, due to the changes above a small change was made into _scanvi.py file. Navigate (SEE ON server if this change is mandatory):
`nano ~/miniconda3/envs/PYenv/lib/python3.10/site-packages/scvi/model/_scanvi.py`

4. Following code, which extracts the scvi-model from the returned tuple, was added to row 189, right after the `from_scvi_model` function:
 - `scvi_model = scvi_model[1]`

