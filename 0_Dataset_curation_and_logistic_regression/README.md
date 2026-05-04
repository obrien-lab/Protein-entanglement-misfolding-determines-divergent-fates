
Use the jupyter notebook `Statistical_association.ipynb` and the data files to reproduce the data curation, logistic regression and protein selection for the [Meadow 2024](https://doi.org/10.1016/j.mcpro.2024.100791) dataset.

Use the jupyter notebook `Statistical_association_other_datasets.ipynb` and the data files to reproduce the data curation and logistic regression for the three additional datasets obtained from [Wagner 2011](https://doi.org/10.1074/mcp.M111.013284) and [Kim 2011](https://doi.org/10.1016/j.molcel.2011.08.025).

The data files include:
| Files | Contents |
| ------ | ------ |
| `AF_Human_shape.xlsx` | Asphericity values for AlphaFold proteins in the dataset, computed by [HullRad](http://52.14.70.9/HullRadV10.1.py) |
| `AF4_knots.dat` | Knotted protein obtained from [AlphaKnot 2.0](https://alphaknot.cent.uw.edu.pl/) |
| `Dataset_S1.xlsx` | Human fibroblast birthdating Ubq-MS data |
| `Human_AF_combined_20250614.csv` | Entanglement information for human proteins obtained from the AlphaFold structures. Details can be found [here](https://www.nature.com/articles/s41467-025-66236-3) |
| `Human_Avg_pLDDTs.csv` | Average pLDDT values for AlphaFold proteins |
| `ipi.HUMAN.xrefs` | Legacy IPI to Uniprot mapping data obtained from [this GitHub repository](https://github.com/sacdallago/IPI_to_UniProt). |
| `Other_datasets.xlsx` | Additional datasets obtained from [Wagner 2011](https://doi.org/10.1074/mcp.M111.013284) and [Kim 2011](https://doi.org/10.1016/j.molcel.2011.08.025). |
| `uniprotkb_human_keywords_membrane_2025_03_19.fasta` | Human proteins with annonation keyword "membrane" obtained from Uniprot |
| `uniprotkb_human_monomer_2025_03_18.fasta` | Human proteins with annonation keyword "monomer" obtained from Uniprot |
| `uniprotkb_human_monomeric_2025_03_18.fasta` | Human proteins with annonation keyword "monomeric" obtained from Uniprot |
| `UP000005640_9606.fasta` | Human proteome sequences obtained from Uniprot |

Note that protein selection is a stochastic process and may yield different results when different random seeds are used. Consequently, the original set of proteins selected in this study for MD simulations cannot be exactly reproduced, as the dataset was updated at a later stage.