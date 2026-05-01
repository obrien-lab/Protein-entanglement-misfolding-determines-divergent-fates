import os, sys, string
import numpy as np
import mdtraj as mdt
import freesasa

def calc_freesasa(traj, frame_idx):
    structure = freesasa.Structure()
    topology = traj.topology
    nan_residx_list = []
    for atom in topology.atoms:
        if np.any(np.isnan(traj.xyz[frame_idx][atom.index])):
            print('Warning: skip nan coordinates atom at frame %d, atom %d'%(frame_idx, atom.index))
            nan_residx_list.append(atom.residue.index)
            continue
        structure.addAtom(atom.name, atom.residue.name, atom.residue.resSeq, string.ascii_uppercase[atom.residue.chain.index], *(traj.xyz[frame_idx][atom.index]*10))
    structure.setRadiiWithClassifier(classifier)
    result = freesasa.calc(structure)
    result = np.array([res_sasa_info.total for chain, residue in result.residueAreas().items() for res_sasa_info in residue.values()])
    nan_residx_list = np.unique(nan_residx_list)
    if len(nan_residx_list) > 0:
        result[nan_residx_list] = np.nan
    return result

traj_dir = './aa_traj_last_100ns/'
num_traj = 50
aa_psf_file = os.popen('ls %s/*.psf'%traj_dir).readlines()[0].strip()

freesasa.setVerbosity(1)
freesasa.Parameters().setNSlices(50)
freesasa.Parameters().setNThreads(1)
classifier = freesasa.Classifier('/storage/home/yuj179/work/software/Jwalk/src/Jwalk/naccess.config.txt')

sasa_list = []
for i in range(num_traj):
    traj_file = '%s/%d_prod_aa.dcd'%(traj_dir, i+1)
    traj = mdt.load(traj_file, top=aa_psf_file)
    # sasa = mdt.shrake_rupley(traj, mode='residue')
    sasa = [calc_freesasa(traj, j) for j in range(traj.n_frames)]
    sasa_list.append(sasa)
    print('Done for traj #%d'%(i+1))

np.save('./analysis/aa_SASA.npy', sasa_list)
