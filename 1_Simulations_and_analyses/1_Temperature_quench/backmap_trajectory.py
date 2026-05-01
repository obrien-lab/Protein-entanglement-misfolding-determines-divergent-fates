import os, sys
import mdtraj as mdt
import parmed as pmd

psf_file = sys.argv[1]
cg_traj_file = sys.argv[2]
start_frame_idx = int(sys.argv[3])
end_frame_idx = int(sys.argv[4])
stride = int(sys.argv[5])
aa_ref_pdb_file = sys.argv[6]
out_dir = sys.argv[7]

run_dir = os.getcwd()

aa_ref_pdb_prefix = aa_ref_pdb_file.split('/')[-1].split('.pdb')[0]
aa_ref_pdb_file_abs = os.path.abspath(aa_ref_pdb_file)

traj_name_prefix = cg_traj_file.split('/')[-1].split('.dcd')[0]

tmp_out_dir = '%s/%s/'%(out_dir, traj_name_prefix)

out_traj_file = '%s/%s_aa.dcd'%(out_dir, traj_name_prefix)

os.system('mkdir -p %s/'%(tmp_out_dir))

cg_traj = mdt.load(cg_traj_file, top=psf_file)

if start_frame_idx < 0:
    start_frame_idx += cg_traj.n_frames
if end_frame_idx < 0:
    end_frame_idx += cg_traj.n_frames

if os.path.exists(out_traj_file):
    out_traj = mdt.load(out_traj_file, top='%s/%s.psf'%(out_dir, aa_ref_pdb_prefix))
    start_frame_idx += out_traj.n_frames * stride
else:
    out_traj = None

cg_traj = cg_traj[start_frame_idx:end_frame_idx+1:stride]

if cg_traj.n_frames > 0:
    cg_traj = cg_traj.center_coordinates()
    cg_traj = cg_traj.superpose(cg_traj)

for i in range(cg_traj.n_frames):
    os.chdir(tmp_out_dir)
    cg_pdb = 'tmp.pdb'
    cg_traj[i].save_pdb(cg_pdb, force_overwrite=True)
    print('Backmapping frame #%d:'%(start_frame_idx + i * stride))
    os.system('backmap.py -i %s -c %s -n 1 -p 1'%(aa_ref_pdb_file_abs, cg_pdb))
    if os.path.exists('./tmp_rebuilt.pdb'):
        if out_traj is None:
            out_traj = mdt.load('tmp_rebuilt.pdb')
        else:
            out_traj = out_traj.join(mdt.load('tmp_rebuilt.pdb'), check_topology=False)
        
        # Save AA psf file if not exists
        if not os.path.exists('../%s.psf'%(aa_ref_pdb_prefix)):
            struct = pmd.load_file('tmp_rebuilt.pdb')
            struct.save('../%s.psf'%(aa_ref_pdb_prefix), overwrite=True)
        
        os.system('rm -f ./tmp_rebuilt.pdb')
        os.chdir(run_dir)
    else:
        print('Failed to backmap. Exit.')
        os.chdir(run_dir)
        sys.exit()
    out_traj.save_dcd(out_traj_file, force_overwrite=True)
    
