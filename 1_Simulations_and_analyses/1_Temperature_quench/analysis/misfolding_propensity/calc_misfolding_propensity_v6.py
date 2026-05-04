import os, sys, multiprocessing, time
import numpy as np
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

frame_segment_list = np.arange(0, 20001, 1000) # ns
G_cutoff = 0.005
n_boot = int(1e5)
n_boot_in_perm = int(1e3)
n_perm = int(1e3)
num_proc = 20

######### FUN ############
def bootstrap(boot_fun, data, n_time, boot_fun_args={}):
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    
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
    combined_data = np.array(list(data_1) + list(data_2), dtype=object)
    t0 = perm_fun(perm_stat_fun, np.arange(len(combined_data)), combined_data, len(data_1))
    
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    start_time = time.time()
    print('start permutation test')
    for i in range(boot_data_1.shape[1]):
        d_1 = boot_data_1[:,i]
        d_2 = boot_data_2[:,i]
        combined_data = np.array(list(d_1) + list(d_2), dtype=object)
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
YU_E_list = np.loadtxt('../../../../new/gen_parameter/YU_E.txt', dtype=str)
NU_NE_list = np.loadtxt('../../../../new/gen_parameter/NU_NE.txt', dtype=str)
NU_E_list = np.loadtxt('../../../../new/gen_parameter/NU_E.txt', dtype=str)

label_list = ['YU_E', 'NU_NE', 'NU_E']

label_postfix_list = []
for i in range(len(frame_segment_list)-1):
    label_postfix_list.append('%d to %d \u03bcs'%(frame_segment_list[i]/1000, frame_segment_list[i+1]/1000))
frame_segment_list = np.ceil(frame_segment_list * 1e3 / 0.015 / 5000 /10) # frames

raw_data = {}
data_df_list = []
perm_data_list = []
for idx, data_list in enumerate([YU_E_list, NU_NE_list, NU_E_list]):
    pd_data = []
    p_nn_prot_list = [[] for ln in label_postfix_list]
    p_ent_prot_list = [[] for ln in label_postfix_list]
    p_nnne_prot_list = [[] for ln in label_postfix_list]
    boot_nn_prot_data = [[] for ln in label_postfix_list]
    boot_ent_prot_data = [[] for ln in label_postfix_list]
    boot_nnne_prot_data = [[] for ln in label_postfix_list]
    for uid in data_list:
        work_dir = '../../' + label_list[idx] + '/' + uid + '/analysis/'
        n_traj = np.nan
        p_ent_list = [[] for ln in label_postfix_list]
        p_nn_list = [[] for ln in label_postfix_list]
        p_nnne_list = [[] for ln in label_postfix_list]
        if os.path.exists(work_dir):
            if os.path.exists(work_dir+'msm_data.npz'):
                npz_data = np.load(work_dir+'msm_data.npz')
                coarse_state_centers = npz_data['coarse_state_centers']
                
                native_idx = len(coarse_state_centers)-1
                ent_idx = np.where(coarse_state_centers[:-1,1] >= G_cutoff)[0]
                nn_nent_idx = np.where(coarse_state_centers[:-1,1] < G_cutoff)[0]
                #native_idx = nent_idx[-1]
                #nn_nent_idx = nent_idx[:-1]
                
                for iframe in range(len(label_postfix_list)):
                    frame_sel = np.arange(frame_segment_list[iframe], frame_segment_list[iframe+1], 1, dtype=int)
                    meta_dtrajs = npz_data['meta_dtrajs'][:, frame_sel]
                
                    n_traj = meta_dtrajs.shape[0]
                
                    for i_traj in range(n_traj):
                        nf_native = len(np.where(meta_dtrajs[i_traj,:] == native_idx)[0])
                        nf_ent = 0
                        for i in ent_idx:
                            nf_ent += len(np.where(meta_dtrajs[i_traj,:] == i)[0])
                        nf_nn_nent = 0
                        for i in nn_nent_idx:
                            nf_nn_nent += len(np.where(meta_dtrajs[i_traj,:] == i)[0])
                        p_nn = 1 - nf_native/len(meta_dtrajs[i_traj,:])
                        p_nn_list[iframe].append(p_nn)
                        p_ent = nf_ent/len(meta_dtrajs[i_traj,:])
                        p_ent_list[iframe].append(p_ent)
                        p_nn_nent = nf_nn_nent/len(meta_dtrajs[i_traj,:])
                        p_nnne_list[iframe].append(p_nn_nent)

        pd_data.append([uid, n_traj])
        raw_data[uid] = [p_nn_list, p_ent_list, p_nnne_list]
        
        for iframe in range(len(label_postfix_list)):
            if len(p_nn_list[iframe]) == 0:
                p_nn_str = ''
                p_nn_prot_list[iframe].append(np.nan)
                boot_nn_prot_data[iframe].append(np.nan*np.ones(n_boot))
            else:
                mean = np.mean(p_nn_list[iframe])
                boot_stat = bootstrap(np.mean, np.array(p_nn_list[iframe]), n_boot)
                ub = np.percentile(boot_stat, 97.5)
                lb = np.percentile(boot_stat, 2.5)
                p_nn_str = '%.2f [%.2f, %.2f]'%(mean, lb, ub)
                p_nn_prot_list[iframe].append(mean)
                boot_nn_prot_data[iframe].append(boot_stat)

            if len(p_ent_list[iframe]) == 0:
                p_ent_str = ''
                p_ent_prot_list[iframe].append(np.nan)
                boot_ent_prot_data[iframe].append(np.nan*np.ones(n_boot))
            else:
                mean = np.mean(p_ent_list[iframe])
                boot_stat = bootstrap(np.mean, np.array(p_ent_list[iframe]), n_boot)
                ub = np.percentile(boot_stat, 97.5)
                lb = np.percentile(boot_stat, 2.5)
                p_ent_str = '%.2f [%.2f, %.2f]'%(mean, lb, ub)
                p_ent_prot_list[iframe].append(mean)
                boot_ent_prot_data[iframe].append(boot_stat)
        
            if len(p_nnne_list[iframe]) == 0:
                p_nnne_str = ''
                p_nnne_prot_list[iframe].append(np.nan)
                boot_nnne_prot_data[iframe].append(np.nan*np.ones(n_boot))
            else:
                mean = np.mean(p_nnne_list[iframe])
                boot_stat = bootstrap(np.mean, np.array(p_nnne_list[iframe]), n_boot)
                ub = np.percentile(boot_stat, 97.5)
                lb = np.percentile(boot_stat, 2.5)
                p_nnne_str = '%.2f [%.2f, %.2f]'%(mean, lb, ub)
                p_nnne_prot_list[iframe].append(mean)
                boot_nnne_prot_data[iframe].append(boot_stat)

            pd_data[-1] += [p_nn_str, p_ent_str, p_nnne_str]

    perm_data_list.append([[np.array(p_nn_prot_list), np.array(boot_nn_prot_data)[:,:,:n_boot_in_perm]],
                           [np.array(p_ent_prot_list), np.array(boot_ent_prot_data)[:,:,:n_boot_in_perm]],
                           [np.array(p_nnne_prot_list), np.array(boot_nnne_prot_data)[:,:,:n_boot_in_perm]]])
    
    uid = 'All'
    d_0 = np.array(pd_data, dtype=object)[:,1].astype(np.float64)
    n_traj = np.sum(d_0[~np.isnan(d_0)].astype(np.int64))
    
    p_nn = np.nanmean(np.array(p_nn_prot_list), axis=1)
    boot_stat = bootstrap_h(np.nanmean, np.array(boot_nn_prot_data).T, boot_fun_args={'axis':0}).T
    ub = np.nanpercentile(boot_stat, 97.5, axis=1)
    lb = np.nanpercentile(boot_stat, 2.5, axis=1)
    p_nn_str = ['%.2f [%.2f, %.2f]'%(p_nn[i], lb[i], ub[i]) for i in range(len(p_nn))]
    
    p_ent = np.nanmean(np.array(p_ent_prot_list), axis=1)
    boot_stat = bootstrap_h(np.nanmean, np.array(boot_ent_prot_data).T, boot_fun_args={'axis':0}).T
    ub = np.percentile(boot_stat, 97.5, axis=1)
    lb = np.percentile(boot_stat, 2.5, axis=1)
    p_ent_str = ['%.2f [%.2f, %.2f]'%(p_ent[i], lb[i], ub[i]) for i in range(len(p_ent))]
    
    p_nnne = np.nanmean(np.array(p_nnne_prot_list), axis=1)
    boot_stat = bootstrap_h(np.nanmean, np.array(boot_nnne_prot_data).T, boot_fun_args={'axis':0}).T
    ub = np.percentile(boot_stat, 97.5, axis=1)
    lb = np.percentile(boot_stat, 2.5, axis=1)
    p_nnne_str = ['%.2f [%.2f, %.2f]'%(p_nnne[i], lb[i], ub[i]) for i in range(len(p_nnne))]
    
    pd_data.append(['All', n_traj])
    columns = ['Uniprot', 'Num trajs']
    for i in range(len(p_nn_str)):
        pd_data[-1] += [p_nn_str[i], p_ent_str[i], p_nnne_str[i]]
        columns += ['Non-native fraction (%s)'%label_postfix_list[i], 
                    'Entangled fraction (%s)'%label_postfix_list[i], 
                    'Non-native non-entangled fraction (%s)'%label_postfix_list[i]]
    
    data_df = pd.DataFrame(pd_data, columns=columns)
    data_df_list.append(data_df)

np.save('raw_data_G%.3f_v6.npy'%(G_cutoff), raw_data)

# Permutation test
p_val_list = []
columns = ['']
for i in range(len(perm_data_list[0][0][0])): # frame segments
    columns += ['Non-native fraction (%s)'%label_postfix_list[i], 
                'Entangled fraction (%s)'%label_postfix_list[i], 
                'Non-native non-entangled fraction (%s)'%label_postfix_list[i]]
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

with pd.ExcelWriter('./misfolding_propensity_G%.3f_v6.xlsx'%(G_cutoff)) as writer:
    for idx, df in enumerate(data_df_list):
        df.to_excel(writer, sheet_name=label_list[idx], index=False)
    pval_df.to_excel(writer, sheet_name='p-value', index=False)


# Plot bars
fig = plt.figure(figsize=(9, 5))
plt.subplots_adjust(left=0.1, wspace=0.3, hspace=0.5, top=0.9, bottom=0.1)
label_list = ['YU & E', 'NU & NE', 'NU & E']
for i_ax in range(9):
    ax = fig.add_subplot(3,3,1+i_ax)
    bar_data = []
    bar_yerr = []
    for j in range(len(label_list)):
        data_list = data_df_list[j].dropna().iloc[:,i_ax-9].to_list()
        mean_list = np.array([float(d.split()[0]) for d in data_list])
        lb_list = np.array([float(d.split()[1][1:-1]) for d in data_list])
        ub_list = np.array([float(d.split()[2][:-1]) for d in data_list])
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
    y_label = data_df_list[0].columns[i_ax-9].split('(')[0].strip()
    if len(y_label) > 25:
        words = y_label.split()
        words[0] += '\n'
        y_label = ' '.join(words)
    ax.set_ylabel(y_label)
    ax.set_ylim([0,1])
    ax.set_title(label_postfix_list[int(i_ax/3)-3])

fig.savefig('./misfolding_propensity_G%.3f_v6.svg'%(G_cutoff))

# plot time series
fig = plt.figure(figsize=(6, 5))
plt.subplots_adjust(left=0.2, wspace=0.4, hspace=0.5, top=0.9, bottom=0.1)
label_list = ['YU & E', 'NU & NE', 'NU & E']
color_list = ['red', 'blue', 'orange']
for i_ax in range(6):
    ax = fig.add_subplot(3,2,1+i_ax)
    if i_ax % 2 == 0:
        for j in range(len(label_list)):
            data_list = data_df_list[j].dropna().iloc[-1,2+int(i_ax/2)::3].to_list()
            mean_list = np.array([float(d.split()[0]) for d in data_list])
            lb_list = np.array([float(d.split()[1][1:-1]) for d in data_list])
            ub_list = np.array([float(d.split()[2][:-1]) for d in data_list])
        
            t_span = np.array([int(l.split()[2]) for l in label_postfix_list])
            ax.plot(t_span, mean_list, c=color_list[j], label=label_list[j])
        
            poly_xy = np.zeros((mean_list.shape[0]*2, 2))
            poly_xy[:mean_list.shape[0], 0] = t_span
            poly_xy[:mean_list.shape[0], 1] = lb_list
            poly_xy[mean_list.shape[0]:mean_list.shape[0]*2, 0] = np.flip(t_span)
            poly_xy[mean_list.shape[0]:mean_list.shape[0]*2, 1] = np.flip(ub_list)
            patch = matplotlib.patches.Polygon(poly_xy, edgecolor=None, 
                                               facecolor=color_list[j], 
                                               alpha=0.2)
            ax.add_patch(patch)
        y_label = data_df_list[0].columns[2+int(i_ax/2)].split('(')[0].strip()
        if len(y_label) > 25:
            words = y_label.split()
            words[0] += '\n'
            y_label = ' '.join(words)
        ax.set_ylabel(y_label)
        ax.set_ylim([0,1])
    else:
        # p-values vs t
        data_list = data_df_list[-1].dropna().iloc[:,1+int(i_ax/2)::3].to_numpy()
        for j in range(len(data_list)):
            label = data_df_list[-1].dropna().iloc[j,0]
            label = label.replace('_', ' & ')
            ax.plot(t_span, data_list[j], c=color_list[j], label=label)
        y_label = 'p-values ' + data_df_list[0].columns[2+int(i_ax/2)].split('(')[0].strip()
        if len(y_label) > 25:
            words = y_label.split()
            words[1] += '\n'
            y_label = ' '.join(words)
        ax.set_ylabel(y_label)
        ax.set_ylim([0,1.5])
        ax.set_yscale('symlog', linthresh=1e-3, base=10)
        ax.axhline(0.05, linestyle='--', lw=1, color='k')

    ax.set_xlabel('Simulation time (\u03bcs)')
    ax.legend()

fig.savefig('./misfolding_propensity_vs_t_G%.3f_v6.svg'%(G_cutoff))
