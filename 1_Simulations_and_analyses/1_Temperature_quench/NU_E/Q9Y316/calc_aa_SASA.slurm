#!/bin/bash
#SBATCH -J SA_Q9Y316
#SBATCH --partition=basic
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 7-00:00:00
#SBATCH --account=epo2_cr_default

cd $SLURM_SUBMIT_DIR

python ../../calc_aa_sasa.py
