| Folder | Contents |
| ------ | -------- |
| `YU_E` | Subfolders containing input files for the 10 young-ubiquitinated and entangled proteins |
| `NU_NE` | Subfolders containing input files for the 10 non-ubiquitinated and non-entangled proteins |
| `NU_E` | Subfolders containing input files for the 10 non-ubiquitinated and entangled proteins |

Within each subfolder, you will find:
- A `setup` folder containing:
  - The AlphaFold structure (`.pdb` format)
  - The domain definition file (`domain_def.dat`)
- A SLURM script for submitting the computational job to a GPU cluster

The code used for these jobs, along with usage instructions, can be found [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/opt_nscal.py).

| File | Contents |
| ------ | -------- |
| `YU_E.txt` | List of Uniprot IDs of the 10 young-ubiquitinated and entangled proteins |
| `NU_NE.txt` | List of Uniprot IDs of the 10 non-ubiquitinated and non-entangled proteins |
| `NU_E.txt` | List of Uniprot IDs of the 10 non-ubiquitinated and entangled proteins |
