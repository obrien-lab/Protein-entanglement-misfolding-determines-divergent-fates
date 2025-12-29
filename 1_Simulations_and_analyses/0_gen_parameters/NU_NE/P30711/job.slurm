#!/bin/bash
#SBATCH -J P30711
#SBATCH --partition=mgc-nih
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --gres=gpu:2
#SBATCH -t 14-00:00:00
#SBATCH --account=epo2_nih

cd $SLURM_SUBMIT_DIR

opt_nscal.py -i setup/P30711.pdb -d setup/domain_def.dat -t 310 -n 0 -p 2 
