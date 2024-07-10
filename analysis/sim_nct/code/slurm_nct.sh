#!/bin/bash
#SBATCH --job-name=scDist_nct
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=philnicol740@gmail.com
#SBATCH --mem=64gb
#SBATCH --output=../logs/scDist_nct.log

module unload R/4.3.2b
module load gcc/9.2.0
module load R/4.1

Rscript num_ct.R
