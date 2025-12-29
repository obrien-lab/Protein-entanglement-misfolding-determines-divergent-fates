#!/usr/bin/bash
#SBATCH -J a_A6NDU8_G
#SBATCH --partition=basic
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 14-00:00:00
#SBATCH --account=epo2_cr_default

cd $SLURM_SUBMIT_DIR

calc_cont_synth_G_vs_T_v2.py 0 20
