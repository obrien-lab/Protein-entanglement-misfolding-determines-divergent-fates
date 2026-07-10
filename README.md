# Protein entanglement misfolding determines divergent fates
This repository contains scripts and input files to reproduce the work described in the manuscript titled **"Protein entanglement misfolding determines divergent fates: proteasomal degradation or persistence in near-native misfolded states."**

The folder `0_Dataset_curation_and_logistic_regression` contains scripts and input files for:
- Processing the human proteome birth-dating dataset  
- Integrating entanglement datasets with the birth-dating dataset  
- Performing logistic regression to test statistical associations while controlling confounding factors such as protein length and oligomer status.
- Selecting proteins for downstream molecular dynamics (MD) simulations  
- Analyzing additional Ubq-MS datasets

This part required 200MB of storage, 1 CPU hour with 16 GB of RAM, no GPUs were required.

The folder `1_Simulations_and_analyses` contains scripts and input files for:
- Running MD simulations  
- Performing analyses on the resulting simulation trajectories  

This part required 2.2TB of storage, 10^5 CPU hours with 128 GB of RAM and 10^4 GPU hours. 

All simulation trajectories are publicly available on [CyVerse](https://data.cyverse.org/dav-anon/iplant/projects/NCEMS/working-groups/protein-misfolding-aging/data/Protein-entanglement-misfolding-determines-divergent-fates/).

The refolding trajectories of the temperature quench simulations for all 30 proteins are also available on MDRepo, with the following IDs:
| YU-E proteins | MDRepo ID |
| -------- | -------- |
| Q0PNE2 | [MDR00022226](https://mdrepo.org/explore/22226) |
| P16152 | [MDR00022184](https://mdrepo.org/explore/22184) |
| P29218 | []() |
| P00491 | [MDR00022173](https://mdrepo.org/explore/22173) |
| P19623 | []() |
| A8MXV4 | [MDR00022171](https://mdrepo.org/explore/22171) |
| P04350 | [MDR00022177](https://mdrepo.org/explore/22177) |
| P31150 | []() |
| O95394 | [MDR00022172](https://mdrepo.org/explore/22172) |
| P52888 | [MDR00022278](https://mdrepo.org/explore/22278) |

| NU-NE proteins | MDRepo ID |
| -------- | -------- |
| Q8WV22 | [MDR00022239](https://mdrepo.org/explore/22239) |
| A6NDU8 | [MDR00022148](https://mdrepo.org/explore/22148) |
| P30711 | [MDR00022223](https://mdrepo.org/explore/22223) |
| Q9BU89 | [MDR00022242](https://mdrepo.org/explore/22242) |
| Q13825 | [MDR00022282](https://mdrepo.org/explore/22282) |
| Q6NVY1 | [MDR00022232](https://mdrepo.org/explore/22232) |
| P16520 | [MDR00022188](https://mdrepo.org/explore/22188) |
| Q8IV38 | [MDR00022235](https://mdrepo.org/explore/22235) |
| P02774 | [MDR00022174](https://mdrepo.org/explore/22174) |
| Q12996 | [MDR00022281](https://mdrepo.org/explore/22281) |

| NU-E proteins | MDRepo ID |
| -------- | -------- |
| P07738 | [MDR00022181](https://mdrepo.org/explore/22181) |
| Q9UIV1 | [MDR00022274](https://mdrepo.org/explore/22274) |
| Q9UBP6 | [MDR00022280](https://mdrepo.org/explore/22280) |
| Q9Y316 | [MDR00022275](https://mdrepo.org/explore/22275) |
| Q63HM1 | [MDR00022276](https://mdrepo.org/explore/22276) |
| P39748 | []() |
| Q6DKJ4 | [MDR00022229](https://mdrepo.org/explore/22229) |
| Q9HB40 | [MDR00022279](https://mdrepo.org/explore/22279) |
| Q96C11 | [MDR00022277](https://mdrepo.org/explore/22277) |
| Q9H6R3 | []() |

Please refer to the `README` files within each folder for detailed instructions and additional information.

### Required python packages:
Ensure that the following Python packages are installed before running the scripts in this repository. For additional scripts from other repositories, please refer to their respective usage pages for any extra package requirements.
```text
biopython==1.81
freesasa==2.2.1
matplotlib==3.5.3
mdtraj==1.9.7
msmtools==1.2.6
networkx==2.6.3
numpy==1.21.6
pandas==1.3.5
ParmEd==3.4.3
pyEMMA==2.5.12
scipy==1.7.3
statsmodels==0.13.5
```

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
