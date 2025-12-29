import numpy as np
import os, sys, time

start_id = 0
total_extend = 50
njob_mgc_max = 15
job_prefix = 'P_P16152_'

####################################

def get_njob(allocation_key, job_prefix):
    word = os.popen('squeue -u yuj179 | grep %s | grep %s | wc -l'%(allocation_key[:9], job_prefix[:8])).readline().strip()
    njob = int(word)
    word = os.popen('squeue -u yuj179 | grep %s | grep " CA " | grep %s | wc -l'%(allocation_key[:9], job_prefix[:8])).readline().strip()
    njob -= int(word)
    return njob

####################################
localtime = '_'.join(time.asctime(time.localtime(time.time())).split())
f = open('Monitor_%s.log'%localtime, 'w')
f.write('Monitor starts at %s\n'%(time.asctime(time.localtime(time.time()))))
f.write('Start id: %d\n'%start_id)
f.close()
while start_id < total_extend:
    njob_mgc = get_njob('mgc-nih', job_prefix)
    print(njob_mgc)
    n_submit = np.min([njob_mgc_max - njob_mgc, total_extend - start_id])
    if n_submit > 0:
        f = open('Monitor_%s.log'%localtime, 'a')
        f.write('Submit %d jobs to mgc-nih at %s\n'%(n_submit, time.asctime(time.localtime(time.time()))))
        lines = os.popen('python extend_sim_roar.py %d %d'%(start_id, n_submit)).readlines()
        for line in lines:
            f.write(line)
        start_id += n_submit
        f.write('Start id: %d\n'%start_id)
        f.close()
    
    time.sleep(2)
    

