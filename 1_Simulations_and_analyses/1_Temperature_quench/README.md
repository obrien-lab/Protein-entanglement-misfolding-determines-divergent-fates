| Folder | Contents |
| ------ | -------- |
| `YU_E` | Subfolders containing input files for the 10 young-ubiquitinated and entangled proteins |
| `NU_NE` | Subfolders containing input files for the 10 non-ubiquitinated and non-entangled proteins |
| `NU_E` | Subfolders containing input files for the 10 non-ubiquitinated and entangled proteins |
| `analysis/misfolding_propensity` | Python script for computing misfolding propensity using the last 500 ns, 300 ns, and 100 ns of the simulation trajectories |
| `analysis/Q_norm` | Python script for computing Q<sub>norm</sub> using the last 500 ns, 300 ns, and 100 ns of the simulation trajectories |

Within each protein subfolder, you will find:

- A `setup` folder containing:
  - AlphaFold structure (`.pdb` format)
  - Secondary structural elements definition file (`secondary_struc_defs.txt`)
  - Domain definition file (`domain_def.dat`)
  - Coarse-grained (CG) model files (`.psf`, `.cor`, `.prm`, `.top`)

- `Tq.cntrl`  
  Simulation setup parameters for temperature-quench simulations. Details can be found [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/temperature_quenching.py).

- `job.slurm`  
  SLURM script for submitting the simulation job to a GPU cluster. The underlying code and usage instructions are available [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/temperature_quenching.py).

- `analysis_Q.slurm`  
  SLURM script for submitting the analysis job for the order parameter **Q** to a CPU cluster. The code and usage instructions are available [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/calc_native_contact_fraction.pl).

- `analysis_G.slurm`  
  SLURM script for submitting the analysis job for the order parameter **G** to a CPU cluster. The code and usage instructions are available [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/calc_entanglement_number_v2.pl).

- `analysis_chirality.slurm`  
  SLURM script for submitting the analysis job for the order parameter **K** to a CPU cluster. The code and usage instructions are available [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/calc_chirality_number.pl).

- An `analysis` folder containing:
  - `get_order_parameters.py`  
    Python script to parse order parameters Q, G, and K from text files and save them as NumPy arrays
  - `get_rm_traj_id.py`  
    Python script to identify trajectories containing mirror-image structures, which are removed in downstream analyses
  - `MSM_sample.py`  
    Python script to cluster simulation structures into metastable states, extract representative structures, and plot −ln(P) and state maps
  - `MSM_sample.cntrl`  
    Configuration file for `MSM_sample.py`. Update file paths as needed when running on your system
  - `get_MSTS.py`  
    Python script to plot state probabilities as a function of simulation time

