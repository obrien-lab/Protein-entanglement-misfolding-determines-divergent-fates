#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import parmed as pmd
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.labelsize'] = 'large'
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 4
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['xtick.labelsize'] = 'medium'
matplotlib.rcParams['ytick.labelsize'] = 'medium'
matplotlib.rcParams['legend.fontsize'] = 'medium'
matplotlib.rcParams['figure.dpi'] = 600

mutant_type = 'A6NDU8'

psf_file = os.popen('ls ../setup/*_ca.psf').readlines()[0].strip()
cor_file = os.popen('ls ../setup/*_ca.cor').readlines()[0].strip()

max_length = len(pmd.load_file(psf_file).atoms)
    
# Get number of native contact
cutoff = 0.8 # in nm
min_interval = 4
struct = pmd.load_file(psf_file)
coor = pmd.load_file(cor_file)
struct.coordinates = coor.coordinates
native_contact = []
for i in range(0, len(struct.atoms)-min_interval):
    coor_1 = np.array([struct[i].xx, struct[i].xy, struct[i].xz])
    for j in range(i+min_interval, len(struct.atoms)):
        coor_2 = np.array([struct[j].xx, struct[j].xy, struct[j].xz])
        dist = np.sum((coor_1-coor_2)**2)**0.5
        if (dist < 10*cutoff):
            native_contact.append([i, j])
num_nc = len(native_contact)
print('# of native contacts: %d'%num_nc)
    
# Get Entanglement trajectory
G_list_0 = list(np.load(mutant_type+'_Entanglement.npy', allow_pickle=True))
G_list = [np.sum(g[:,:5], axis=1)/num_nc for g in G_list_0]
    
# NC length trajectories
T_list_0 = list(np.load(mutant_type+'_T.npy', allow_pickle=True))
NCL_list = [t[:,1] for t in T_list_0]
stage_list = [t[:,2] for t in T_list_0]

G_list_end = []
for traj_idx in range(len(G_list)):
    idx_list = np.where(stage_list[traj_idx] == 3)[0]
    dL_list = NCL_list[traj_idx][idx_list][1:] - NCL_list[traj_idx][idx_list][:-1]
    idx_list_1 = np.append(np.where(dL_list == 1)[0], len(idx_list)-1)
    G_list_end.append(G_list[traj_idx][idx_list][idx_list_1])
    
    #idx_list = np.where(stage_list[traj_idx] == 5)[0][-1]
    G_list_end[-1] = np.append(G_list_end[-1], G_list[traj_idx][-1])

# plot
fig_width = 5
fig_height = fig_width*0.3
fig = plt.figure(figsize=(fig_width, fig_height))
plt.subplots_adjust(top=0.9, bottom=0.3, left=0.2, right=0.9)

ax = fig.add_subplot(1,1,1)
for G_list_1 in G_list_end:
    ax.plot(np.arange(1,max_length+1,dtype=int), G_list_1, '-', color=np.array([135, 206, 250])/255, lw=0.5)
ax.plot(np.arange(1,max_length+1,dtype=int), np.mean(G_list_end, axis=0), '-b', lw=2.0, label='Average')
ax.set_xlabel('Nascent chain length (AA)')
ax.set_ylabel('$\mathsf{G}$')
ax.legend()

ax.set_ylim((-0.0010683760978281499, 0.022435898054391147))

fig.savefig('G_vs_length.svg')

data = np.vstack((np.arange(1,max_length+1), np.array(G_list_end)))
data = np.vstack((data, np.mean(G_list_end, axis=0)))
data = data.T
df = pd.DataFrame(data, columns=['Nascent chain length (AA)']+['Trajectory #%d'%(i+1) for i in range(len(G_list_end))]+['Average'])

df.to_excel('G_vs_length.xlsx', index=False)
