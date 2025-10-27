#!/bin/bash
#SBATCH --job-name=snakemake_controller
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/snakemake_controller_%j.out
#SBATCH --partition=high
#SBATCH --account=ACCOUNTNAME

#run this scripts in running directory (usually one directory above /scripts)

source ~/.bashrc
conda activate snakemake_env

# low partition recommended for snakemake command. 
# run this snakemake controller script on high partitiona with #SBATCH so that the controller doesn't get kicked off. 
# if the controller stops, the whole Snakefile stops. 
snakemake -s scripts/Snakefile_step1_3_FINAL.py \
    --executor slurm \
    --jobs 20 \
    --use-conda \
    --rerun-incomplete \
    --printshellcmds \
    --default-resources mem_mb=4096 runtime=600 slurm_partition=low