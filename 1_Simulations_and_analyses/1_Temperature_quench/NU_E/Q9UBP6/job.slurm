#!/bin/bash
#SBATCH -J Tq_Q9UBP6
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --gres=gpu:1
#SBATCH -t 14-00:00:00
#SBATCH --account=epo2_cr_default

cd $SLURM_SUBMIT_DIR

temperature_quenching.py -f Tq.cntrl
