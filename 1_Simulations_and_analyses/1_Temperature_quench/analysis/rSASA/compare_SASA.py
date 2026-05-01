import os, sys, multiprocessing, time
import numpy as np
import parmed as pmd
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.labelsize'] = 'small'
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 4
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['xtick.labelsize'] = 'x-small'
matplotlib.rcParams['ytick.labelsize'] = 'x-small'
matplotlib.rcParams['legend.fontsize'] = 'x-small'
matplotlib.rcParams['figure.dpi'] = 600

last_nframes_list = np.array([-100]) # ns
G_cutoff = 0.005
n_boot = int(1e5)
n_boot_in_perm = int(1e3)
n_perm = int(1e3)
num_proc = 8
res_sel_mask = ':ILE,VAL,LEU,PHE,CYS,MET,ALA,GLY,TRP'

######### FUN ############
def bootstrap(boot_fun, data, n_time, boot_fun_args={}):
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    
    if len(data) == 0:
        return np.nan * np.ones(n_time)
    
    idx_list = np.arange(len(data))
    
    print('start bootsrapping')
    for i in range(n_time):
        sample_idx_list = np.random.choice(idx_list, len(idx_list))
        new_data = data[sample_idx_list]
        pool_list.append(pool.apply_async(boot_fun, args=(new_data, ), kwds=boot_fun_args))
    pool.close()
    pool.join()
    boot_stat = np.array([p.get() for p in pool_list])
    return boot_stat
    
def bootstrap_h(boot_fun, data, boot_fun_args={}):
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    
    idx_list = np.arange(data.shape[1])
    
    print('start bootsrapping')
    for i in range(data.shape[0]):
        sample_idx_list = np.random.choice(idx_list, len(idx_list))
        new_data = data[i, sample_idx_list]
        pool_list.append(pool.apply_async(boot_fun, args=(new_data, ), kwds=boot_fun_args))
    pool.close()
    pool.join()
    boot_stat = np.array([p.get() for p in pool_list])
    return boot_stat
    
def permutation_test(perm_stat_fun, data_1, data_2, num_perm, perm_fun):
    combined_data = np.array(list(data_1) + list(data_2))
    t0 = perm_fun(perm_stat_fun, np.arange(len(combined_data)), combined_data, len(data_1))
    
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    start_time = time.time()
    print('start permutation test')
    for i in range(num_perm):
        perm_idx_list = np.random.permutation(np.arange(len(combined_data)))
        pool_list.append(pool.apply_async(perm_fun, (perm_stat_fun, perm_idx_list, combined_data, len(data_1))))
    pool.close()
    pool.join()
    t_dist = [p.get() for p in pool_list]
    p = 0
    for t in t_dist:
        if t >= t0:
            p += 1
    p = (p+1)/(num_perm+1)
    used_time = time.time() - start_time
    print('%.2fs'%used_time)
    return p
    
def permutation_test_h(perm_stat_fun, data_1, data_2, boot_data_1, boot_data_2, num_perm, perm_fun):
    combined_data = np.array(list(data_1) + list(data_2))
    t0 = perm_fun(perm_stat_fun, np.arange(len(combined_data)), combined_data, len(data_1))
    
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    start_time = time.time()
    print('start permutation test')
    for i in range(boot_data_1.shape[1]):
        d_1 = boot_data_1[:,i]
        d_2 = boot_data_2[:,i]
        combined_data = np.array(list(d_1) + list(d_2))
        for j in range(num_perm):
            perm_idx_list = np.random.permutation(np.arange(len(combined_data)))
            pool_list.append(pool.apply_async(perm_fun, (perm_stat_fun, perm_idx_list, combined_data, len(d_1))))
    pool.close()
    pool.join()
    t_dist = np.array([p.get() for p in pool_list])
    t_dist[np.isnan(t_dist)] = []
    p = (len(np.where(t_dist >= t0)[0])+1) / (len(t_dist)+1)
    used_time = time.time() - start_time
    print('%.2fs'%used_time)
    return p
    
def perm_fun(perm_stat_fun, perm_idx_list, combined_data, length_1):
    d_1 = perm_stat_fun(combined_data[perm_idx_list[:length_1]])
    d_2 = perm_stat_fun(combined_data[perm_idx_list[length_1:]])
    return np.abs(d_1-d_2)
    
########## MAIN #############
YU_E_list = np.loadtxt('../../../gen_parameter/YU_E.txt', dtype=str)
NU_E_list = np.loadtxt('../../../gen_parameter/NU_E.txt', dtype=str)

label_list = ['YU_E', 'NU_E']

label_postfix_list = []
for i in range(len(last_nframes_list)):
    if i == len(last_nframes_list)-1:
        label_postfix_list.append('last %d ns'%(-last_nframes_list[i]))
    else:
        label_postfix_list.append('last %d to last %d ns'%(-last_nframes_list[i], -last_nframes_list[i+1]))
last_nframes_list = np.ceil(last_nframes_list * 1e3 / 0.015 / 5000) # frames

raw_data = {}
data_df_list = []
perm_data_list = []
for idx, data_list in enumerate([YU_E_list, NU_E_list]):
    pd_data = []
    p_ent_prot_list = [[] for ln in last_nframes_list]
    p_nent_prot_list = [[] for ln in last_nframes_list]
    boot_ent_prot_data = [[] for ln in last_nframes_list]
    boot_nent_prot_data = [[] for ln in last_nframes_list]
    for uid in data_list:
        work_dir = '../../' + label_list[idx] + '/' + uid + '/analysis/'
        n_traj = np.nan
        SASA_ent_list = [[] for ln in last_nframes_list]
        SASA_nent_list = [[] for ln in last_nframes_list]
        SASA_native_list = [[] for ln in last_nframes_list]
        if os.path.exists(work_dir):
            if os.path.exists(work_dir+'msm_data.npz'):
                npz_data = np.load(work_dir+'msm_data.npz')
                coarse_state_centers = npz_data['coarse_state_centers']

                aa_pdb_file = os.popen('ls '+work_dir+'../setup/*.pdb').readlines()[0].strip()
                aa_struct = pmd.load_file(aa_pdb_file)
                cg_psf_file = os.popen('ls '+work_dir+'../setup/*.psf').readlines()[0].strip()
                cg_struct = pmd.load_file(cg_psf_file)
                for i in range(len(cg_struct.residues)):
                    cg_struct.residues[i].name = aa_struct.residues[i].name
                    cg_struct.residues[i].atoms[0].name = 'CA'
                    cg_struct.residues[i].atoms[0].element = 6
                res_sel_idx = list(pmd.amber.AmberMask(cg_struct, res_sel_mask).Selected())

                n_traj = npz_data['meta_dtrajs'].shape[0]
                
                skip_traj_list = []
                file_object = open(work_dir+'MSM_sample.cntrl','r')
                for line in file_object:
                    if line.startswith('skip_traj_list'):
                        words = line.split('=')
                        skip_traj_list = [int(w) for w in words[1].strip().split()]
                mask = np.ones(len(skip_traj_list)+n_traj, dtype=bool)
                mask[skip_traj_list] = False
                
                # Get SASA trajectory
                SASA_list = np.load(work_dir+'aa_SASA.npy', allow_pickle=True)[mask]
                SASA_list = np.sum(SASA_list[:,:,res_sel_idx], axis=-1)
                
                native_idx = len(coarse_state_centers)-1
                ent_idx = np.where(coarse_state_centers[:-1,1] >= G_cutoff)[0]
                nent_idx = np.where(coarse_state_centers[:-1,1] < G_cutoff)[0]
                
                for iframe in range(len(last_nframes_list)):
                    if iframe == len(last_nframes_list)-1:
                        frame_sel = np.arange(last_nframes_list[iframe], 0, 1, dtype=int)
                    else:
                        frame_sel = np.arange(last_nframes_list[iframe], last_nframes_list[iframe+1], 1, dtype=int)
                    meta_dtrajs = npz_data['meta_dtrajs'][:, frame_sel]
                
                    SASA_list_0 = SASA_list[:, frame_sel]
                    
                    native_SASA = []
                    ent_SASA = []
                    nent_SASA = []
                    for it, md in enumerate(meta_dtrajs):
                        idx_list = np.where(md == native_idx)[0]
                        if len(idx_list) > 0:
                            native_SASA.append(np.nanmean(SASA_list_0[it, idx_list]))
                        idx_list = np.where(np.isin(md, ent_idx))[0]
                        if len(idx_list) > 0:
                            ent_SASA.append(np.nanmean(SASA_list_0[it, idx_list]))
                        idx_list = np.where(np.isin(md, nent_idx))[0]
                        if len(idx_list) > 0:
                            nent_SASA.append(np.nanmean(SASA_list_0[it, idx_list]))
                    
                    SASA_ent_list[iframe] = ent_SASA
                    SASA_nent_list[iframe] = nent_SASA
                    SASA_native_list[iframe] = native_SASA                    

        pd_data.append([uid, n_traj])
        raw_data[uid] = [[SASA_native_list, SASA_ent_list, SASA_nent_list]]
        # bootstrap
        boot_native_data = []
        boot_ent_data = []
        boot_nent_data = []
        for iframe in range(len(last_nframes_list)):
            boot_native_data.append(bootstrap(np.mean, np.array(SASA_native_list[iframe]), n_boot))
            boot_ent_data.append(bootstrap(np.mean, np.array(SASA_ent_list[iframe]), n_boot))
            boot_nent_data.append(bootstrap(np.mean, np.array(SASA_nent_list[iframe]), n_boot))
        raw_data[uid].append([boot_native_data, boot_ent_data, boot_nent_data])
        
        # compute mean and 95% CI
        for iframe in range(len(last_nframes_list)):
            if len(SASA_native_list[iframe]) == 0:
                p_ent_str = ''
                p_ent_prot_list[iframe].append(np.nan)
                boot_ent_prot_data[iframe].append(np.nan*np.ones(n_boot))
                p_nent_str = ''
                p_nent_prot_list[iframe].append(np.nan)
                boot_nent_prot_data[iframe].append(np.nan*np.ones(n_boot))
            else:
                if len(SASA_ent_list[iframe]) == 0:
                    p_ent_str = ''
                    p_ent_prot_list[iframe].append(np.nan)
                    boot_ent_prot_data[iframe].append(np.nan*np.ones(n_boot))
                else:
                    ratio = np.mean(SASA_ent_list[iframe]) / np.mean(SASA_native_list[iframe])
                    boot_stat_ratio = boot_ent_data[iframe] / boot_native_data[iframe]
                    ub = np.percentile(boot_stat_ratio, 97.5)
                    lb = np.percentile(boot_stat_ratio, 2.5)
                    p_ent_str = '%.2f [%.2f, %.2f]'%(ratio, lb, ub)
                    p_ent_prot_list[iframe].append(ratio)
                    boot_ent_prot_data[iframe].append(boot_stat_ratio)
                if len(SASA_nent_list[iframe]) == 0:
                    p_nent_str = ''
                    p_nent_prot_list[iframe].append(np.nan)
                    boot_nent_prot_data[iframe].append(np.nan*np.ones(n_boot))
                else:
                    ratio = np.mean(SASA_nent_list[iframe]) / np.mean(SASA_native_list[iframe])
                    boot_stat_ratio = boot_nent_data[iframe] / boot_native_data[iframe]
                    ub = np.percentile(boot_stat_ratio, 97.5)
                    lb = np.percentile(boot_stat_ratio, 2.5)
                    p_nent_str = '%.2f [%.2f, %.2f]'%(ratio, lb, ub)
                    p_nent_prot_list[iframe].append(ratio)
                    boot_nent_prot_data[iframe].append(boot_stat_ratio)
                
            pd_data[-1] += [p_ent_str, p_nent_str]

    perm_data_list.append([[np.array(p_ent_prot_list), np.array(boot_ent_prot_data)[:,:,:n_boot_in_perm]],
                           [np.array(p_nent_prot_list), np.array(boot_nent_prot_data)[:,:,:n_boot_in_perm]]])
    
    uid = 'All'
    d_0 = np.array(pd_data, dtype=object)[:,1].astype(np.float64)
    n_traj = np.sum(d_0[~np.isnan(d_0)].astype(np.int64))
    
    p_ent = np.nanmean(np.array(p_ent_prot_list), axis=1)
    boot_stat = bootstrap_h(np.nanmean, np.array(boot_ent_prot_data).T, boot_fun_args={'axis':0}).T
    ub = np.nanpercentile(boot_stat, 97.5, axis=1)
    lb = np.nanpercentile(boot_stat, 2.5, axis=1)
    p_ent_str = ['%.2f [%.2f, %.2f]'%(p_ent[i], lb[i], ub[i]) for i in range(len(p_ent))]
    
    p_nent = np.nanmean(np.array(p_nent_prot_list), axis=1)
    boot_stat = bootstrap_h(np.nanmean, np.array(boot_nent_prot_data).T, boot_fun_args={'axis':0}).T
    ub = np.nanpercentile(boot_stat, 97.5, axis=1)
    lb = np.nanpercentile(boot_stat, 2.5, axis=1)
    p_nent_str = ['%.2f [%.2f, %.2f]'%(p_nent[i], lb[i], ub[i]) for i in range(len(p_nent))]

    pd_data.append(['All', n_traj])
    columns = ['Uniprot', 'Num trajs']
    for i in range(len(p_ent_str)):
        pd_data[-1] += [p_ent_str[i], p_nent_str[i]]
        columns += ['Entangled fraction (%s)'%label_postfix_list[i],  
                    'Non-entangled fraction (%s)'%label_postfix_list[i]]
    
    data_df = pd.DataFrame(pd_data, columns=columns)
    data_df_list.append(data_df)

np.save('SASA_ratio_raw_data_G%.3f_v4.npy'%(G_cutoff), raw_data)

# Permutation test
p_val_list = []
columns = ['']
for i in range(len(perm_data_list[0][0][0])): # frame segments
    columns += ['Entangled fraction (%s)'%label_postfix_list[i], 
                'Non-entangled fraction (%s)'%label_postfix_list[i]]
for pi in range(len(perm_data_list)-1):
    for pj in range(pi+1, len(perm_data_list)):
        p_val_list.append([label_list[pi] + ' vs. ' + label_list[pj]])
        for i in range(len(perm_data_list[0][0][0])): # frame segments
            for j in range(len(perm_data_list[0])): # states
                data_1 = perm_data_list[pi][j][0][i]
                boot_data_1 = perm_data_list[pi][j][1][i]
                data_2 = perm_data_list[pj][j][0][i]
                boot_data_2 = perm_data_list[pj][j][1][i]
                p_val = permutation_test_h(np.nanmean, data_1, data_2, boot_data_1, boot_data_2, n_perm, perm_fun)
                p_val_list[-1].append(p_val)

pval_df = pd.DataFrame(p_val_list, columns=columns)

with pd.ExcelWriter('./SASA_ratio_G%.3f_v4.xlsx'%(G_cutoff)) as writer:
    for idx, df in enumerate(data_df_list):
        df.to_excel(writer, sheet_name=label_list[idx], index=False)
    pval_df.to_excel(writer, sheet_name='p-value', index=False)
        
# Plot
fig = plt.figure(figsize=(6, 1.7))
plt.subplots_adjust(left=0.1, wspace=0.3, hspace=0.5, top=0.9, bottom=0.1)
label_list = ['YU & E', 'NU & E']
for i_ax in range(2):
    ax = fig.add_subplot(1,2,1+i_ax)
    bar_data = []
    bar_yerr = []
    for j in range(len(label_list)):
        data_list = data_df_list[j].dropna().iloc[:,2+i_ax].to_list()
        mean_list = np.array([float(d.split()[0]) if d != '' else np.nan for d in data_list])
        lb_list = np.array([float(d.split()[1][1:-1]) if d != '' else np.nan for d in data_list])
        ub_list = np.array([float(d.split()[2][:-1]) if d != '' else np.nan for d in data_list])
        y_err = np.array([mean_list-lb_list, ub_list-mean_list])
        ax.errorbar(np.linspace(j-0.4, j+0.4, len(mean_list[:-1])), mean_list[:-1], 
                    yerr=y_err[:,:-1], fmt='o', markerfacecolor='none', markeredgecolor='gray',
                    markersize=5, ecolor='gray', elinewidth=1, capsize=3)
        bar_data.append(mean_list[-1])
        bar_yerr.append(y_err[:,-1])
    bar_yerr = np.array(bar_yerr).T
    ax.bar(np.arange(len(label_list)), bar_data, width=0.8, bottom=0, yerr=bar_yerr,
           facecolor='grey', edgecolor='black', linewidth=0, tick_label=label_list,
           error_kw = {'ecolor':'black', 'capsize':20, 'elinewidth':2, 'capthick':2},
           alpha=0.2, zorder=10)
    y_label = data_df_list[0].columns[2+i_ax].split('(')[0].strip()
    if len(y_label) > 25:
        words = y_label.split()
        words[0] += '\n'
        y_label = ' '.join(words)
    ax.set_ylabel(y_label)
    ax.set_title(label_postfix_list[int(i_ax/2)])

fig.savefig('./SASA_ratio_G%.3f_v4.svg'%(G_cutoff))
