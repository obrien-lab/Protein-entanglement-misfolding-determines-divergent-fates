import sys, os
import numpy as np

prefix = sys.argv[1]

Q_threshold = 0.4
# Q_threshold = 0.2
K_threshold = 0.5
# K_threshold = 0.6

K_list = np.load('./%s_K.npy'%prefix, allow_pickle=True)
Q_list = np.load('./%s_QBB.npy'%prefix, allow_pickle=True)
Q_list = np.array([Q[:K_list.shape[1],:] for Q in Q_list])
n = Q_list.shape[-1]-1
n_domain = int((np.sqrt(1+8*n)-1)/2)
Q_list = Q_list[:,:,:n_domain]
K_list = K_list[:,:,:n_domain]

Q_mean_list = np.mean(Q_list[:,-134:,:], axis=1) # Last 10 ns
K_mean_list = np.mean(K_list[:,-134:,:], axis=1)

Q_tag = Q_mean_list >= Q_threshold
K_tag = K_mean_list < K_threshold

QK_tag = Q_tag & K_tag
Tag = np.any(QK_tag, axis=-1)

rm_traj_idx_list = np.where(Tag)[0]

for ri in rm_traj_idx_list:
    idx_list = np.array(np.where(QK_tag[ri,:])[0])+1
    print('Traj #%d Domain: '%(ri+1) + str(idx_list))

print(rm_traj_idx_list)


