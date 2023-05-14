#!/bin/bash

### slurm job options: line must begin with '#SBATCH'

#SBATCH --job-name=data    # job name
#SBATCH --mail-type=END,FAIL    # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=erno.hanninen@sund.ku.dk    # email address to receive the notification
#SBATCH -c 5    # number of requested cores
#SBATCH --mem=30gb    # total requested RAM
#SBATCH --time=0-05:00:00               # max. running time of the job, format in D-HH:MM:SS
#SBATCH --output=any_name_%j.log   # Standard output and error log, '%j' gives the job ID

### write your own scripts below

module load java/11.0.15 nextflow

nextflow run data_processing_wf.nf
