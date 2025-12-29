import os, sys

work_dir = '/storage/home/yuj179/mygroup/protein_Ubq/new/Tq/'
folder_list = ['YU_E', 'NU_NE', 'NU_E']

for folder_0 in folder_list:
    original_dir = work_dir+folder_0
    target_dir = './'+folder_0
    folder_list_1 = [d for d in os.listdir(original_dir) if os.path.isdir(os.path.join(original_dir, d))]
    for folder_2 in folder_list_1:
        original_dir_1 = os.path.join(original_dir, folder_2)
        target_dir_1 = os.path.join(target_dir, folder_2)
        os.system('mkdir -p %s'%target_dir_1)
        os.system('cp -r %s/setup %s/'%(original_dir_1, target_dir_1))
        os.system('cp -r %s/analysis %s/'%(original_dir_1, target_dir_1))
        os.system('cp %s/job.slurm %s/'%(original_dir_1, target_dir_1))
        os.system('cp %s/Tq.cntrl %s/'%(original_dir_1, target_dir_1))
        os.system('cp %s/analysis_chirality.slurm %s/'%(original_dir_1, target_dir_1))
        os.system('cp %s/analysis_G.slurm %s/'%(original_dir_1, target_dir_1))
        os.system('cp %s/analysis_Q.slurm %s/'%(original_dir_1, target_dir_1))
