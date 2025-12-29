#!/bin/bash
#SBATCH -J Tq_P16520
#SBATCH --partition=mgc-nih
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --gres=gpu:1
#SBATCH -t 10-00:00:00
#SBATCH --account=epo2_nih

cd $SLURM_SUBMIT_DIR

temperature_quenching.py -f Tq.cntrl
