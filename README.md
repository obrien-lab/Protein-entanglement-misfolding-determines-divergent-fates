# Protein entanglement misfolding determines divergent fates
This repository contains scripts and input files to reproduce the work described in the manuscript titled **"Protein entanglement misfolding determines divergent fates: proteasomal degradation or persistence in near-native misfolded states."**

The folder `0_Dataset_curation_and_logistic_regression` contains scripts and input files for:
- Processing the human proteome birth-dating dataset  
- Integrating entanglement datasets with the birth-dating dataset  
- Performing logistic regression to test statistical associations  
- Selecting proteins for downstream molecular dynamics (MD) simulations  

This part required 200MB of storage, 1 CPU hour with 16 GB of RAM, no GPUs were required.

The folder `1_Simulations_and_analyses` contains scripts and input files for:
- Running MD simulations  
- Performing analyses on the resulting simulation trajectories  

This part required 2.2TB of storage, 10^5 CPU hours with 128 GB of RAM and 10^4 GPU hours.

Please refer to the `README` files within each folder for detailed instructions and additional information.

### Contact: Prof. Ed O'Brien epo2@psu.edu

### How to cite: 
<details>
<summary>Click to expand BibTeX</summary>

```bibtex
@article {Jiang2026.04.15.718748,
    author = {Jiang, Yang and Jain, Anushka and Ghaemmaghami, Sina and O{\textquoteright}Brien, Edward P.},
    title = {Protein entanglement misfolding determines divergent fates: proteasomal degradation or persistence in near-native misfolded states},
    elocation-id = {2026.04.15.718748},
    year = {2026},
    doi = {10.64898/2026.04.15.718748},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2026/04/16/2026.04.15.718748},
    eprint = {https://www.biorxiv.org/content/early/2026/04/16/2026.04.15.718748.full.pdf},
    journal = {bioRxiv}
}
```
</details>

This work was supported by the National Science Foundation National Synthesis Center for the Emergence of Molecular and Cellular Sciences NCEMS (DBI-2335029)