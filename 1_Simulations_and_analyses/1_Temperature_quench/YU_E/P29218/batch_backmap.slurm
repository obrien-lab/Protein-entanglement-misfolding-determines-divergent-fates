#!/usr/bin/env python3
#SBATCH -J bm_P29218
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

def analysis(traj):
    try:
        os.system("python ../../backmap_trajectory.py %s %s %d %d %d %s %s"%(psf_file, traj, start_frame_idx, end_frame_idx, stride, aa_ref_pdb_file, out_dir))
    except Exception as e:
        traceback.print_exc()

nproc = 20
num_traj = 50
start_frame_idx = -1333 # last 100 ns
end_frame_idx = -1 
stride = 1 
aa_ref_pdb_file = os.popen('ls setup/*.pdb').readlines()[0].strip()
os.chdir(os.environ['SLURM_SUBMIT_DIR'])

psf_file = os.popen('ls ./setup/*.psf').readlines()[0].strip()

out_dir = './aa_traj_last_100ns/'

traj_list = [str(i) for i in range(1,num_traj+1)]

os.system("mkdir -p %s"%out_dir)

pool = multiprocessing.Pool(nproc)
for traj_dir in traj_list:
    pool.apply_async(analysis, ('./traj/'+traj_dir+'_prod.dcd',))

while True:
    time.sleep(10)
    if len(pool._cache) == 0:
        print('All Done.')
        break

pool.close()
pool.join()
