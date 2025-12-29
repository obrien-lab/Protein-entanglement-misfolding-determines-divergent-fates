import numpy as np
import os, sys, random

walltime = '10-00:00:00'
n_rep = 1
n_traj = 50
sim_time = 1e-6 # microsec
submit_start_id = int(sys.argv[1])
submit_num = int(sys.argv[2])

psffile = 'prot_l294.psf'
ncrstfile = 'prot_l294_dissociation_final.ncrst'
fffile = '/storage/home/yuj179/mygroup/protein_Ubq/post_translation/A6NDU8/1-50/setup/A6NDU8_clean_nscal1_fnn1_go_bt.xml'
corfile = '/storage/home/yuj179/mygroup/protein_Ubq/post_translation/A6NDU8/1-50/setup/A6NDU8_clean_ca.cor'
sec_struct_def_file = '/storage/home/yuj179/mygroup/protein_Ubq/post_translation/A6NDU8/1-50/setup/secondary_struc_defs.txt'
alpha = 1
timestep = 0.015
sim_steps = int(np.ceil(sim_time*1e9/alpha*1e3 / timestep / 5000) * 5000)

for ei in range(submit_start_id, submit_start_id+submit_num):
    idx_traj = int(ei/n_rep)+1
    idx_rep = int(ei-(idx_traj-1)*n_rep)+1
    outfile = '%d/traj_%d_%d.out'%(idx_traj, idx_traj, idx_rep)
    tag_finish = False
    if os.path.exists(outfile):
        f = open(outfile)
        lines = f.readlines()
        f.close()
        last_line = lines[-1].strip()
        if last_line.startswith('Done'):
            tag_finish = True
        elif not last_line.startswith('Time') and os.path.getsize(outfile) != 0:
            end_step = int(last_line.split()[1])
            if end_step >= sim_steps:
                tag_finish = True
        
    if not tag_finish:
        os.chdir('%d/'%(idx_traj))
        
        rand = int(random.random()*1e7)
        
        f = open('job_%d.slurm'%(idx_rep), 'w')
        f.write('''#!/bin/bash
#SBATCH -J P_A6NDU8_'''+str(idx_traj)+'''_'''+str(idx_rep)+'''
#SBATCH --partition=mgc-nih
#SBATCH --mail-type=END
#SBATCH -o %j.o # Name of stdout output file
#SBATCH -e %j.e # Name of stderr error file
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --gres=gpu:1
#SBATCH -t '''+walltime+'''
#SBATCH --account=epo2_nih

cd $SLURM_SUBMIT_DIR

post_trans_single_run.py '''+psffile+''' '''+ncrstfile+''' '''+fffile+''' 310 1 traj_'''+str(idx_traj)+'''_'''+str(idx_rep)+''' '''+str(rand)+''' '''+str(sim_steps)+''' '''+sec_struct_def_file+''' 1.1 '''+corfile+''' 0
''')
        f.close()
        
        os.system('sbatch job_%d.slurm'%(idx_rep))
        print('Submit %d %d'%(idx_traj, idx_rep))
        
        os.chdir('../')
    else:
        print('Skip %d %d finished job'%(idx_traj, idx_rep))
