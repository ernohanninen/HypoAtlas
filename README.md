1. Create Python environment:
`conda create --name PYenv`

2. Activate:
`conda activate PYenv`

3. Install packages thru conda:
`conda install pymde==0.1.18 scvi-tools=0.19.0 jupyter==1.0.0 nbconvert==7.2.5`

4. Install packages thru pip:
`pip install nbparameterise==0.5 scib==1.0.4`

5. Install s4vectors:
`conda install -c bioconda bioconductor-s4vectors==0.36.0`

6. Install singlecellexperiment:
`conda install -c bioconda bioconductor-singlecellexperiment==1.20.0`

7. Install seurat
`conda install -c bioconda r-seurat==4.3.0`

8. Install scgen, trvae, travaep:
`pip install scgen trvae trvaep`


Or one liner: 
`conda create --name PYenv && conda activate PYenv && conda install pip pymde==0.1.18 scvi-tools=0.19.0 jupyter==1.0.0 nbconvert==7.2.5 && pip install nbparameterise==0.5 scib==1.0.4 && conda install -c bioconda bioconductor-s4vectors==0.36.0 && conda install -c bioconda bioconductor-singlecellexperiment==1.20.0 && conda install -c bioconda r-seurat==4.3.0 && pip install scgen trvae trvaep`


The conda environments used to benchmark different data integration tools:
- ComBat: PYenv, for ComBat algorithm scib uses scanpy implementation. Therefore using PYenv containing SCANPY package to benchmark ComBat


- mnnpy: mmn_env, created by running:
`conda create --name mnn_env jupyter==1.0.0 mnnpy scanpy pip && conda activate mnn_env && pip install scib && conda deactivate`


- desc: desc_env, created by running:
`conda create --name desc_env pip && conda activate desc_env && pip install jupyter scib desc tensorflow && conda deactivate`


- saucie: saucie_env, created by running:
1. installing jupyter and scib, and cloning the repo:
`conda create --name saucie_env pip && conda activate saucie_env && pip install jupyter scib && git clone https://github.com/KrishnaswamyLab/SAUCIE`
2. Due to package conflicts, the requirements.txt file was edited:
`nano SAUCIE/requirements.txt`
3. The requirements.txt file used:
`tensorflow
fcsparser
fcswrite
numpy
pandas
scikit-learn
matplotlib
`
4. Install the dependencies:
`pip install -r SAUCIE/requirements.txt && conda deactivate`


- bbknn: bbknn_env, created by running:
`conda create --name bbknn_env jupyter==1.0.0 pip && conda activate bbknn_env && pip install scib scib[bbknn] && conda deactivate`


- harmony: harmony_env, created by running:
`conda create --name harmony_env jupyter==1.0.0 pip && conda activate harmony_env && pip install scib harmony-pytorch && conda deactivate`

- scanorama: scanorama_env, created by running:
`conda create --name scanorama_env python-annoy pip && conda activate scanorama_env && pip install scib jupyter && git clone https://github.com/brianhie/scanorama.git && cd scanorama/ && python setup.py install --user && conda deactivate`

- SCVI: PYenv, scvi-tools are already installed in PYenv, therefore using it to run scvi-algorithm

- SCANVI: PYenv, scvi-tools are already installed in PYenv, therefore using it to run scanvi-algorithm

- scGen: PYenv, scgen is installed in PYenv, therefore using it to run scGen algorithm

- trvaep: scgen is installed in PYenv, therefore using it to run scGen algorithm

- trVAE: the traVAE installation to PYenv didn't succeed, to avoid conflicts in PYenv an external environment for traVAE was created by running:
`conda create --name trvae_env jupyter==1.0.0 pip && conda activate trvae_env && pip install scib trvae && conda deactivate`















_________________________________________________________________
conda scanvi
conda create --name scanvi_env jupyter==1.0.0 pip && conda activate scanvi_env && pip install scib scib[scanvi] && conda deactivate


__________________________________________________________
1. Create env for scanorama:
`conda create --name scanorama_env jupyter==1.0.0 pip && conda activate scanorama_env && pip install scib[scanorama] && conda deactivate`

AttributeError: module 'collections' has no attribute 'MutableSet'






_________________________________________________________
1.  Check harmony










Have nbparametrize in the python env PYenv

pip install nbparameterise
add scvi-tools to pyenv




trash:
1. Create CONDA environment:
`conda create --name integration`

2. Activate it: 
`conda activate integration`

3. Install SCVI-tools
`conda install -c conda-forge scvi-tools=0.19.0`

4. Install pymde:
`conda install -c pytorch -c conda-forge pymde`


