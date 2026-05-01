| Folder | Contents |
| ------ | -------- |
| `YU_E` | Subfolders containing input files for the 10 young-ubiquitinated and entangled proteins |
| `NU_NE` | Subfolders containing input files for the 10 non-ubiquitinated and non-entangled proteins |
| `NU_E` | Subfolders containing input files for the 10 non-ubiquitinated and entangled proteins |

Within each subfolder, you will find:
- A `setup` folder containing:
  - The AlphaFold structure (`.pdb` format) downloaded from AlphaFold database.
  - The domain definition file (`domain_def.dat`) created based on the CATH/Pfam domains.
- A SLURM script for submitting the computational job to a GPU cluster

The code used for these jobs, along with usage instructions, can be found [here](https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/opt_nscal.py).

| File | Contents |
| ------ | -------- |
| `YU_E.txt` | List of Uniprot IDs of the 10 young-ubiquitinated and entangled proteins |
| `NU_NE.txt` | List of Uniprot IDs of the 10 non-ubiquitinated and non-entangled proteins |
| `NU_E.txt` | List of Uniprot IDs of the 10 non-ubiquitinated and entangled proteins |

These lists of proteins were obtained from running [`Statistical_association.ipynb`](../../0_Dataset_curation_and_logistic_regression/Statistical_association.ipynb).

Go to each folder and the subfolders, use `sbatch job.slurm` to submit the computational job to HPC. Please adjust the resource configuration in job.slurm accordingly. These jobs will find the minimal nscal value that stablize the protein native structure and generate the psf, cor, prm and secondary structure elements. These files will be used in [`1_Temperature_quench`](../1_Temperature_quench/)