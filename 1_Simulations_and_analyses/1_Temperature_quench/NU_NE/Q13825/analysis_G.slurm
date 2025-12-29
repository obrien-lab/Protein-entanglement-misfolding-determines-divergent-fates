#!/usr/bin/env python3
#SBATCH -J G_Q13825
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

def analysis(cor_file, traj, out, start_frame_id):
	os.system("calc_entanglement_number_v2.pl -i %s -t %s -o %s -b %d -c 0.6 -r 1"%(cor_file, traj, out, start_frame_id))

np = 20
num_traj = 50
start_frame_idx = 0
os.chdir(os.environ['SLURM_SUBMIT_DIR'])

psf_file = os.popen('ls setup/*.psf').readlines()[0].strip()
cor_file = os.popen('ls setup/*.cor').readlines()[0].strip()

traj_list = [str(i) for i in range(1,num_traj+1)]

os.system("mkdir -p analysis/G")

pool = multiprocessing.Pool(np)
for traj_dir in traj_list:
    if start_frame_idx < 0:
        num_frame = mdt.load('traj/'+traj_dir+'_prod.dcd', top=psf_file).n_frames
        start_frame_id = num_frame + start_frame_idx +1
    else:
        start_frame_id = start_frame_idx +1
    pool.apply_async(analysis, (cor_file, 'traj/'+traj_dir+'_prod.dcd', 'analysis/G', start_frame_id))

while True:
    time.sleep(10)
    if len(pool._cache) == 0:
        print('All Done.')
        break

pool.close()
pool.join()
