#!/bin/bash
#SBATCH -J job_moni_A6NDU8
#SBATCH --partition=mgc-nih
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 14-00:00:00
#SBATCH --account=epo2_nih

cd $SLURM_SUBMIT_DIR

python job_monitor_roar.py
