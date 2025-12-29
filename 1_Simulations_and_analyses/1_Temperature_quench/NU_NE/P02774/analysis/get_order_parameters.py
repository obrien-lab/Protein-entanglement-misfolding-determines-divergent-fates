#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import parmed as pmd
import mdtraj as mdt

################################# Arguments ###################################
n_traj = 50
mutant_type = sys.argv[1]
if len(sys.argv) == 2:
    frame_offset = 0
elif len(sys.argv) == 3:
    frame_offset = int(sys.argv[2])

po_dir = '../'

cor_file = os.popen('ls %s/setup/*_ca.cor'%po_dir).readlines()[0].strip()
psf_file = os.popen('ls %s/setup/*_ca.psf'%po_dir).readlines()[0].strip()

sec_def = os.popen('ls %s/setup/secondary_*'%po_dir).readlines()[0].strip()
dom_def = os.popen('ls %s/setup/domain_*'%po_dir).readlines()[0].strip()

################################# Functions ###################################    
def get_Qbb():
    print('QBB:')
    qbb_list = [];
    for i in range(n_traj):
        fo = open(po_dir+'/analysis/qbb/qbb_'+str(i+1)+'_prod.dat', 'r')
        C = fo.readlines()
        fo.close()
        C = [C[k].strip().split() for k in range(1,len(C))]
        C = np.array(C, dtype=np.float32)
        # c_start = C.shape[1]-1
        c_start = 0
        c_end = C.shape[1]
        Q_ts = C[:,c_start:c_end].reshape((len(C),c_end-c_start))
        qbb_list.append(Q_ts[frame_offset:])
        print('  Traj %d Done'%(i+1))
    np.save('%s_QBB_%d.npy'%(mutant_type, frame_offset), qbb_list)
    
def get_entanglement():
    global n_traj, rep_per_traj, co_dir, max_length, cor_file, po_dir, mutant_type
    print('G:')
    G_number_list = [];
    for i in range(n_traj):
        G_list_1 = []
        f = open(po_dir+'/analysis/G/G_'+str(i+1)+'_prod.dat', 'r')
        C = f.readlines()
        f.close()
        C1 = [C[k].strip().split() for k in range(8,len(C))]
        C1 = np.array(C1, dtype=np.float32)
        #print(C1.shape)
        G = C1[frame_offset:,:]
        G_number_list.append(G)
        print('  Traj %d Done %d'%(i+1, G_number_list[-1].shape[0]))   
    np.save('%s_G_%d.npy'%(mutant_type, frame_offset), G_number_list)
    
def get_K():
    print('K:')
    K_list = [];
    for i in range(n_traj):
        fo = open(po_dir+'/analysis/chirality/K_'+str(i+1)+'_prod.dat', 'r')
        C = fo.readlines()
        fo.close()
        C = [C[k].strip().split() for k in range(1,len(C))]
        C = np.array(C, dtype=np.float32)
        # c_start = C.shape[1]-1
        c_start = 0
        c_end = C.shape[1]
        K_ts = C[:,c_start:c_end].reshape((len(C),c_end-c_start))
        K_list.append(K_ts[frame_offset:])
        print('  Traj %d Done'%(i+1))
    np.save('%s_K_%d.npy'%(mutant_type, frame_offset), K_list)
    
################################## MAIN #######################################
structure = pmd.load_file(psf_file)
max_length = len(structure.residues)

# Get Qbb trajectory
get_Qbb()

# Get Entanglement trajectory
get_entanglement()

# Get Chirality trajectory
get_K()
