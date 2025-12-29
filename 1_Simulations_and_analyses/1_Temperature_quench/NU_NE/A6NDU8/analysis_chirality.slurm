#!/usr/bin/env python3
#SBATCH -J K_A6NDU8
#SBATCH --partition=basic
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH --account=epo2_cr_default

import getopt, os, time, multiprocessing, random, math, traceback, io, sys
import mdtraj as mdt

def analysis(cor_file, dom_def, sec_def, traj, out_dir):
    try:
        os.system("calc_chirality_number.pl -i "+cor_file+" -d "+dom_def+" -s "+sec_def+" -t "+traj+" -o "+out_dir)
    except Exception as e:
        traceback.print_exc()

nproc = 20
num_traj = 50

os.chdir(os.environ['SLURM_SUBMIT_DIR'])

psf_file = os.popen('ls setup/*.psf').readlines()[0].strip()
cor_file = os.popen('ls setup/*.cor').readlines()[0].strip()
sec_def = os.popen('ls setup/secondary_*').readlines()[0].strip()
dom_def = os.popen('ls setup/domain_def*').readlines()[0].strip()

traj_list = [str(i) for i in range(1,num_traj+1)]

os.system("mkdir -p analysis/chirality")

pool = multiprocessing.Pool(nproc)
for traj_dir in traj_list:
    pool.apply_async(analysis, (cor_file, dom_def, sec_def, 'traj/'+traj_dir+'_prod.dcd', 'analysis/chirality/'))

while True:
    time.sleep(10)
    if len(pool._cache) == 0:
        print('All Done.')
        break

pool.close()
pool.join()
