
Use the jupyter notebook `Statistical_association.ipynb` and the data files to reproduce the data curation, logistic regression and protein selection.

The data files include:
| Files | Contents |
| ------ | ------ |
| `AF_Human_shape.xlsx` | Asphericity values for AlphaFold proteins in the dataset, computed by [HullRad](http://52.14.70.9/HullRadV10.1.py) |
| `AF4_knots.dat` | Knotted protein obtained from [AlphaKnot 2.0](https://alphaknot.cent.uw.edu.pl/) |
| `Avg_pLDDTs.csv` | Average pLDDT values for AlphaFold proteins |
| `Dataset_S1.xlsx` | Human fibroblast birthdating Ubq-MS data |
| `Human_AF_combined_20250614.csv` | Entanglement information for human proteins obtained from the AlphaFold structures. Details can be found [here](https://www.nature.com/articles/s41467-025-66236-3) |
| `uniprotkb_human_keywords_membrane_2025_03_19.fasta` | Human membrane proteins obtained from Uniprot |
| `UP000005640_9606.fasta` | Human protein sequences obtained from Uniprot |

Note that protein selection is a stochastic process and may yield different results when different random seeds are used. Consequently, the original set of proteins selected in this study for MD simulations cannot be exactly reproduced, as the dataset was updated at a later stage.