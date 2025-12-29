#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects
import msmtools

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['xtick.labelsize'] = 'large'
matplotlib.rcParams['ytick.labelsize'] = 'large'
matplotlib.rcParams['legend.fontsize'] = 'large'

color_map_hex_list = ['#001219', '#005f73', '#0a9396', '#94d2bd', '#e9d8a6', '#ee9b00', '#ca6702', '#bb3e03', '#6d0004']
float_list = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]

################################# Arguments ###################################
end_t = 2000 # in ns
dt = 0.015/1000
nsave = 5000
dt = dt*nsave # in ns
n_traj = 50
n_window = 200
mutant_type_list = ['Tq']
mutant_label_list = ['Temperature quench']
start_idx = 1
end_idx = 1
num_points_plot = 1000
n_boot = 10000
num_proc = 4
if_extend_dtraj = True
if_boot = True

skip_traj = [[]]

################################# Functions ###################################

def boot_fun(data, n_states, job_id):
    PPT = np.zeros((data.shape[1], n_states))
    for md in data:
        for i in range(data.shape[1]):
            PPT[i,md[i]]+=1
    PPT /= len(data)
    # print('Bootstrapping done for %d'%(job_id+1))
    return PPT

def bootstrap(boot_fun, data, n_boot):
    global num_proc, n_states
    
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    
    idx_list = np.arange(len(data))
    
    start_time = time.time()
    print('start bootsrapping')
    for i in range(n_boot):
        sample_idx_list = np.random.choice(idx_list, len(idx_list))
        new_data = data[sample_idx_list,:]
        pool_list.append(pool.apply_async(boot_fun, (new_data, n_states, i, )))
    pool.close()
    pool.join()
    boot_stat = [p.get() for p in pool_list]
    used_time = time.time() - start_time
    print('%.2fs'%used_time)
    return boot_stat

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]
    
def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmap = matplotlib.colors.LinearSegmentedColormap('my_cmap', segmentdata=cdict, N=256)
    return cmap

################################## MAIN #######################################
color_map = get_continuous_cmap(color_map_hex_list, float_list=float_list)
meta_dtrajs = np.load('./msm_data.npz', allow_pickle=True)['meta_dtrajs']

max_T_len = int(np.ceil(end_t/dt))
print('Max trajectory length: %d'%max_T_len)
interval = int(max_T_len / num_points_plot)
sample_idx = [max_T_len-1-i*interval for i in range(int(max_T_len/interval), -1, -1)]
if sample_idx[0] != 0:
    sample_idx = [0] + sample_idx

n_states = 0
for md in meta_dtrajs:
    if n_states < np.max(md):
        n_states = np.max(md)
n_states += 1   

if if_extend_dtraj:
    k = 0
    meta_dtrajs_extended = []
    for i, md in enumerate(meta_dtrajs):
        (N, be) = np.histogram(md[-n_window:], bins=np.arange(-0.5, n_states, 1))
        meta_dtraj_last = np.argwhere(N == np.max(N))[0][0]
        if md[-1] != n_states-1 and len(md) < max_T_len:
            print(i, md[-1], meta_dtraj_last, len(md))
            k+=1
            if i < n_traj * end_idx:
                skip_traj[0].append(i)
            else:
                skip_traj[1].append(i)
        else:
            mde = []
            for j in range(max_T_len):
                if j >= len(md): 
                    if md[-1] == n_states-1:
                        state_0 = n_states-1
                    else:
                        state_0 = meta_dtraj_last
                else:
                    state_0 = md[j]
                mde.append(state_0)
            meta_dtrajs_extended.append(mde)
    print(k)

    meta_dtrajs_extended = np.array(meta_dtrajs_extended)
else:
    all_data = np.load('MSTS_data.npz', allow_pickle=True)
    meta_dtrajs_extended = all_data['meta_dtrajs_extended']
    skip_traj = all_data['skip_traj']

fig = plt.figure(figsize=(6,8))
plt.subplots_adjust(wspace=0.3, hspace=0.3)
ax_list = []
ylim_list = []
si = 0
MSTS_list = []
boot_stat_list = []
for i_ax, mutant_type in enumerate(mutant_type_list):
    ei = si+n_traj*(end_idx-start_idx+1)-len(skip_traj[i_ax])
    meta_dtrajs = meta_dtrajs_extended[np.arange(si, ei)]
    meta_dtrajs = meta_dtrajs[:,sample_idx]
    si = ei
    
    # MSTS
    PPT = np.zeros((len(sample_idx), n_states))
    t_span = (np.array(sample_idx)+1)*dt
    for md in meta_dtrajs:
        for i in range(len(sample_idx)):
            PPT[i,md[i]]+=1
    PPT /= len(meta_dtrajs)
    MSTS_list.append(PPT)
    
    # bootstrap
    if if_boot:
        boot_stat = bootstrap(boot_fun, meta_dtrajs, n_boot)
    else:
        all_data = np.load('MSTS_data.npz', allow_pickle=True)
        boot_stat = all_data['boot_stat_list'][i_ax]
    boot_stat_list.append(boot_stat)
    PPT_ub = np.zeros(PPT.shape)
    PPT_lb = np.zeros(PPT.shape)
    for i in range(n_states):
        bs_state = []
        for bs in boot_stat:
            bs_state.append(bs[:,i])
        PPT_ub[:,i] = np.percentile(bs_state, 97.5, axis=0)
        PPT_lb[:,i] = np.percentile(bs_state, 2.5, axis=0)
    
    ax = fig.add_subplot(2,1,i_ax+1)
    cmap = color_map
    legend_str = []
    MS_data = []
    columns = []
    for i in range(n_states):
        ax.plot(t_span, PPT[:,i], c=np.array(cmap(i/(n_states-1))))
        legend_str.append('Metastable state #%d'%(i+1))
        y_ub = PPT_ub[:,i]
        y_lb = PPT_lb[:,i]
        poly_xy = np.zeros((PPT.shape[0]*2, 2))
        poly_xy[:PPT.shape[0], 0] = t_span
        poly_xy[:PPT.shape[0], 1] = y_lb
        poly_xy[PPT.shape[0]:PPT.shape[0]*2, 0] = np.flip(t_span)
        poly_xy[PPT.shape[0]:PPT.shape[0]*2, 1] = np.flip(y_ub)
        patch = matplotlib.patches.Polygon(poly_xy, edgecolor=None, 
                                           facecolor=cmap(i/(n_states-1)), 
                                           alpha=0.2)
        ax.add_patch(patch)
            
        MS_data.append(PPT[:,i])
        columns.append('P%d'%(i+1))
        MS_data.append(y_lb)
        columns.append('lower bound (95%CI)')
        MS_data.append(y_ub)
        columns.append('upper bound (95%CI)')
    # Save raw data
    MS_data = np.array(MS_data).T
    df = pd.DataFrame(MS_data, columns=columns, index=t_span)
    df.to_csv('MSTS_%s.csv'%mutant_type)
    #ax.legend(legend_str, loc='best', fontsize='x-small')
    ax.set_xlabel('Simulation Time (ns)')
    ax.set_ylabel('State Probability')
    ax.set_title(mutant_label_list[i_ax])
    ax_list.append(ax)
    ylim_list.append(list(ax.get_ylim()))

ylim_list = np.array(ylim_list)
for i_ax, mutant_type in enumerate(mutant_type_list):
    ax = ax_list[i_ax]
    ax.set_ylim([np.min(ylim_list[:,0]), np.max(ylim_list[:,1])])

fig.savefig('MSTS.svg')

np.savez('MSTS_data.npz', 
         meta_dtrajs_extended = meta_dtrajs_extended,
         skip_traj = skip_traj,
         t_span = t_span,
         MSTS_list = MSTS_list,
         boot_stat_list = boot_stat_list)
