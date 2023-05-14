# Master project

## Description of the project
The purpose of this master's degree project was to use publicly available single-cell and spatial transcriptomics resources to analyze scRNA-seq data generated from our in-house human hypothalamic differentiation protocols. To analyze the cell type composition of the hypothalamic differentiation protocols, we constructed a developing human hypothalamic reference atlas and applied it to predict the cell types of our in-house data. Both available developing human hypothalamus datasets [Herb](https://www.biorxiv.org/content/10.1101/2021.07.20.453090v2) and [Zhou](https://www.sciencedirect.com/science/article/pii/S1934590921004574?via%3Dihub) constructed the reference herein. To ensure the usage of the best available integration method we performed a scRNA-seq data integration benchmarking study. This was achieved using a custom-made benchmarking pipeline. Results from the benchmarking study suggest scANVI was the best-performing method. The scANVI training process was optimized by fine-tuning the hyperparameters and adjusting the cell type annotations. The scANVI reference model was used to predict cell types from our in-house data. To validate our in-house scRNA-seq data from hypothalamic differentiation protocols it was spatially aligned towards a spatial dataset from the human fetal neural tube. 

This readme file contains a detailed description of the workflow used to achieve the result of the study. First covering the pre-requests to reproduce the analysis results. Followed by a description of the analysis steps, which in here is divided into four parts: data processing, integration benchmarking, scANVI optimization, and spatial mapping. 

## Pre-requests

### Computing environment
All the commands and analysis presented in this README file were executed on [DanGPU](https://sgn102.pages.ku.dk/a-not-long-tour-of-dangpu/#jupyter-and-rstudio-on-dangpu) HPC-environment, hosted by reNEW, at Copenhagen University. DanGPU uses a slurm workload manager. Some commands / slurm files used during the project might be specific to DanGPU.

### Installations
Clone the repository:

From now on it's assumed the user is located in the ~/HypoAtlas folder

To reproduce the results installations of Conda and Nextflow are required. Here Conda (4.10.1) and Nextflow (22.10.6.5843) were used since they were preinstalled in the computing environment used. 

All the dependencies of the pipelines and scripts are managed using Conda environments. Due to conflicting packages, multiple environments need to be set up. To achieve this run:
```
conda env create -f CondaEnvironments/bbknn_env.yml | conda env create -f CondaEnvironments/bonefight_env.yml | conda env create -f CondaEnvironments/combat_env.yml | conda env create -f CondaEnvironments/desc_env.yml | conda env create -f CondaEnvironments/harmony_env.yml | conda env create -f CondaEnvironments/PYenv.yml | conda env create -f CondaEnvironments/RConvertEnv.yml | conda env create -f CondaEnvironments/Renv.yml | conda env create -f CondaEnvironments/scanorama_env.yml 
```

The two Nextflow pipelines used in this study assume the conda environments are stored in the `~/.conda/envs` folder. If this is not the case the following changes are needed:
1. Adjust the envDir variable accordingly to DataProcessing_pipeline/data_processing_wf.nf file:
```
nano DataProcessing_pipeline/data_processing_wf.nf
```
2. Do the same change to the envDir variable to Benchmarking_pipeline/Scripts/processes.nf file:
```
nano Benchmarking_pipeline/Scripts/processes.nf 
```

### Data

The analysis steps using Nextflow pipelines load the data automatically to the correct location. Since some of the data is not publicly available additional manual work is required. 

scRNA-seq data:
 - The human developing hypothalamic reference atlas presented composes of two datasets [Herb](https://www.biorxiv.org/content/10.1101/2021.07.20.453090v2) and [Zhou](https://www.sciencedirect.com/science/article/pii/S1934590921004574?via%3Dihub). The processed and integrated dataset presented in Zhou et al., paper is not publicly available. Here, a corresponding in-house version of this dataset was used. Before running the analysis this dataset should be placed in `DataProcessing_pipeline/Data/RawZhou` folder. Otherwise, the datasets needed are publicly available. The data processing pipeline loads these datasets from public databases.
 - The reference atlas was used to predict the cell types of our day 16 human hypothalamic differentiation protocol scRNA-seq data. This data, which is not publicly available, should be placed in `Scanvi_notebooks/Data` folder. 
 - In the spatial mapping our in-house hypothalamic differentiation protocol data was supplemented with data from the MiSTR atlas (developing human neural tube). This data, which is not publicly available, should be placed in `Spatial_mapping/Data` folder.

Spatial transcriptomics data:
 - The spatial dataset used in this study is presented [here](https://www.biorxiv.org/content/10.1101/2022.10.24.513487v1.full). Load the dataset using this [link](https://storage.googleapis.com/linnarsson-lab-human/EEL_HE_5week/LBEXP20211113_EEL_HE_5w_970um_RNA_transformed_assigned.parquet) and place it to `Spatial_mapping/Data` folder.


## Analysis workflow

### Processing the raw data
The raw data was quality-filtered and annotated using the Data_processing Nextflow pipeline. The Nextflow pipeline automatically activates the conda environment used to manage the dependencies of the Nextflow process.  The pipeline loads the raw data presented in Zhou et al. (2022) study, and the raw and integrated data presented in Herb et al. (2022) study, from data repositories. The processed and integrated data presented by Zhou et al. (2022) is not publicly available. Therefore before executing the pipeline the in-house annotated dataset (zhou_annotations.RDS), was placed in `DataProcessing_pipeline/Data/RawZhou` folder.

 
To reproduce the results submit the slurm job by running:
```
sbatch DataProcessing_pipeline/data_slurm_job.sh
```

The pipeline outputs a .h5ad file where data from both publications are merged and Jupyter Notebooks are. The output is stored in DataProcessing_pipeline/Data/Processed folder. The merged_zhou_herb.h5ad file, containing the Zhou and Herb data, is the one used for integration benchmarking and when constructing the actual reference atlas, with the best-performing integration method. The executed Notebooks are stored in DataProcessing_pipeline/Output folder.



### scRNA-seq data integration benchmarking study
To get an unbiased view of the performance of different integration methods we developed and used a Nextflow-based benchmarking pipeline. The pipeline uses the output from the data processing pipeline. The benchmarking pipeline is configured to reproduce the results presented in the project report. The pipeline manages the dependencies automatically.

To reproduce the benchmarking results submit the slurm job by running: 
```
sbatch Benchmarking_pipeline/benchmarking_slurm.sh
```

In `Benchmarking_pipeline/Output/integration` folder the pipeline stores a Jupyter Notebook, integration metrics, and integrated data in `.h5ad` or `.rds` format. The Jupyter Notebook contains an overview of each benchmarked integration method. The method output is visualized using a PAGA graph and force-directed graph drawing plot. Based on the benchmarking pipeline scANVI, was the best-performing method.


### Constructing reference atlas with scANVI
The benchmarking study suggests scANVI is the best-performing method. Since it's a deep learning-based method the performance of scANVI is dependent on how precisely its parameters are fine-tuned. However, the scANVI documentation suggests initializing scANVI with a pre-trained scVI-model, so in fact, the parameters are fine-tuned for scVI. In the benchmarking study, we observed cell types causing inconsistent clustering in all integration outputs, and that the cell names between the two publications were not harmonized properly. Therefore scANVI training process was optimized by addressing the cell type issues observed in the benchmarking study and fine-tuning the parameters for scVI. scANVI uses cell type information to guide the integration process. Meaning the cell type adjustments should be done after running scVI, but before running scANVI. The hypothalamic nuclei were annotated only in the Herb dataset, therefore before running scANVI we identified and annotated most of them also from the Zhou dataset. After constructing the reference atlas using scANVI, the reference was used to predict cell types from our in-house day 16 data sequenced from the human hypothalamus differentiation protocol. 

This analysis was done using Jupyter Notebooks since it enables interactive computing, which is essential for scRNA-seq analysis. To run the Jupyter Notebook web interface on DanGPU run the following steps:

1. Reserve resources:
```
srun --gres=gpu:1 --cpus-per-task=18 --mem=300GB --time=0-12:00:00 --pty bash
```

2. PYenv conda environment was used for the following analysis, activate it:
```
conda activate PYenv
```

3. Launch Jupyter:
```
jupyter notebook --no-browser --ip=127.0.0.1 --port=8890
```

4. Let Jupyter run in the current terminal and open a new terminal window

5. From the new terminal window, listen to the port where Jupyter is running :
```
ssh -N -f -L localhost:8890:localhost:8890 bns631@dangpu01fl.unicph.domain
```

7. Copy the URL displayed in the previous terminal and paste it into the browser

Now when we have the Jupyter server running navigate to the `Scanvi_notebooks` folder. This folder contains all notebooks required for constructing the reference dataset and reference-based cell type prediction. Jupyter Notebooks can be executed by opening the desired Notebook and from `Run` select `Run All Cells`. The notebooks are configured to reproduce the result presented in the project report. However, several of the tools used contain randomness meaning that reproducing the exact result might be impossible. Nonetheless, the notebooks should be executed in the following order:
1. `prepare_data_scvi.ipynb`
    - Takes the unintegrated data generated in the data processing pipeline as input
    - Removes NE, Dividing, Blood, and Fibroblasts from the merged dataset
    - Harmonizes the cell type naming between the publications
    - Writes the cell type adjusted AnnData object to the file
2. `scvi_parameter_autotune.ipynb`
    - Reads the output of prepare_data_scvi.ipynb notebook
    - Generates random hyperparameter combinations from the defined search space and evaluates their performance with our data
    - The parameters are ranked based on validation loss
3. `run_scvi.ipynb`
    - Reads the output of prepare_data_scvi.ipynb notebook
    - Uses the best-performing parameter combination identified by the autotune module
    - Integrates data using scVI algorithm
    - Annotates the Tanycytes and Radial glia cell population, which based on marker genes were miss annotated in the Zhou dataset
    - Save the scVI integrated data and scVI model to file
4. `explore_hypothalamic_nuclei.ipynb`
    - Reads the scVI integrated data
    - Hypothalamic nuclei weren't annotated in the Zhou dataset, therefore identify those from scVI integrated data
    - Store the cell type adjusted scVI integrated data to file
5. `run_scanvi.ipynb`
    - Read the scVI integrated data and pre-trained scVI model
    - Adjust the integration using scANVI algorithm, which, unlike scVI, uses available cell type annotations to guide the integration
    - Evaluate scANVI integration using PAGA graph, force-directed graph drawing plot, and marker gene expression
    - Evaluate scANVI model by applying it to predict its cell types
    - Store scANVI integrated data and scANVI model to file
6. d16_prediction.ipynb
    - Read scANVI model
    - Predict cell types from our in-house hypothalamic differentiation protocol data sequenced at differentiation day 16
    - Evaluate prediction results by exploring marker gene expression

Additionally, all these Notebooks were used to create plots for the project report.


### Spatial mapping of day 16 differentiation protocol data
To address how closely our hypothalamic differentiation protocol recapitulates hypothalamic development in vivo, we spatially aligned the data from the protocol to the post-conception week 5 human neural tube spatial dataset. 

The analysis was performed using Jupyter Notebook. Since a different conda environment was used the Jupyter session needs to be re-launched.
1. Reserve resources:
```
srun --gres=gpu:1 --cpus-per-task=18 --mem=300GB --time=0-12:00:00 --pty bash
```

2. bonefight_env conda environment was used for the following analysis, activate it:
```
conda activate bonefight_env
```

3. Launch Jupyter:
```
jupyter notebook --no-browser --ip=127.0.0.1 --port=8890
```

4. Let Jupyter run in the current terminal and open a new terminal window

5. From the new terminal window, listen to the port where Jupyter is running :
```
ssh -N -f -L localhost:8890:localhost:8890 bns631@dangpu01fl.unicph.domain
```

7. Copy the URL displayed in the previous terminal and paste it into the browser


Now when we have the Jupyter server running navigate to the `Spatial_mapping` folder. This folder contains `Hypo_d16_MiSTR_d14_d21_bonefight.ipynb` notebook used to align the scRNA-seq data spatially towards the human neural tube. Open the file and Hypo_d16_MiSTR_d14_d21_bonefight. To execute the notebook from `Run` select `Run All Cells`. The notebook outputs plots used in the report. 


## Limitations
- The benchmarking pipeline contains an implementation of methods which due to installation issues were not benchmarked in this project. To avoid errors, it is required that in 'pipeline.config' file the following variables have the value 'false': 'benchmark_mnn', 'benchmark_saucie', and 'benchmark_conos'. In the future, the code snippets for these algorithms will be removed.
- We observed strange integration results for DESC and trVAE algorithms, indicating either that the tools are implemented incorrectly, or that these are inadequate for complex integration task
- Both the data processing and benchmarking pipeline utilizes scib.pp.read_seurat function to read RDS files to Python. After running the analysis we have encountered errors when using the function. This is most likely since the PYenv (the main environment for this project) has been updated with new packages frequently, and the dependencies for this function are not up-to-date. In case the pipeline is used in future projects, a workaround for this is required, since with the current configuration is it not executable. 










