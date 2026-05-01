| Folder | Contents |
| ------ | -------- |
| `0_gen_parameters` | Input files used to generate coarse-grained (CG) models and force-field parameters for the selected proteins |
| `1_Temperature_quench` | Input files used to run temperature-quench protein refolding simulations for the selected proteins |
| `2_continuous_synthesis` | Input files and analysis scripts used to perform cotranslational protein folding simulations on the ribosome |
| `3_post_translation` | Input files and analysis scripts used to perform post-translational protein folding simulations off the ribosome |


The workflow starts from `0_gen_parameters`, followed by `1_Temperature_quench`, `2_continuous_synthesis` and `3_post_translation`.