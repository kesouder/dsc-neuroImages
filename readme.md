# DSC Capstone Quarter 1 Project: Reproducing Evolutionary Expansion Map vs. Target Maps of Various Brain Structures Using the Neuromaps Toolbox

### Co-Authors: Kevin Huang, Kevin Souder

This repository reproduces analyses from the **[neuromaps](https://www.nature.com/articles/s41592-022-01625-w)** framework — a standardized system for comparing human brain maps across modalities, spaces, and scales.  
Our capstone project uses **neuromaps** to evaluate how surfaced-based **human evolutionary cortical expansion** relates to other cortical features such as cerebral blood volume, myelination, and functional gradients, using the **spin test** for non-parametric spatial null modeling.

## Purpose of this Repository

This repo provides reproducible code for:

- **Fetching** publicly available cortical maps from the `neuromaps.datasets` API (e.g. Hill et al. 2010, Raichle et al. 2010).
- **Resampling** the source map to surface space (fsaverage or civet) and density (< 164k) of target map if necessary.
- **Comparing** the **evolutionary expansion map** (`hill2010:evoexp`) with a set of **target maps** (e.g. `raichle:cbv`) using the **spin permutation test** (Alexander-Bloch et al. 2018).
- **Computing** correlation coefficients (Pearson’s \( r \)) and spin-based \( p \)-values to evaluate statistical significance while accounting for spatial autocorrelation.
- **Plotting** the box plots of distribution of spin test correlations along with observed correlation for each pair of comparison.

## About Neuromaps

**neuromaps** (Markello & Mišić, _Nature Methods_, 2022) is an open-source Python package that standardizes the comparison of brain maps. It provides:

- **Public access** to over 40 published brain maps (e.g., myelin, metabolism, cortical thickness).
- **Standardized coordinate frameworks**, allowing easy transformation across different coordinate systems (fsLR, MNI152, CIVET).
- **Spatial statistics**, various methods for generating spatial nulls for both surface-based and volumemetric-based maps
- **Easy comparison and manipulation**, direct plotting, array-handling, and testing functions applicable for all brain maps.

**Dependencies** For python >= 3.8+, it's possible to directly use `pip` to download neuromaps using `pip install neuromaps`. But it's recommended to install it from
the source directly in order to fetch the latest version and **other required packages**.

```
git clone https://github.com/netneurolab/neuromaps
cd neuromaps
pip install .
```

For plotting, we used `matplotlib`, you can install it in Python using the command

```
pip install matplotlib
```

**IMPORTANT** In order to use tranformations and many other crucial functionalities in `neuromaps`, you must have **[Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)** ready in your path, you can check this by the command

```
wb_command -version
```

**Data Location**
All datasets are automatically fetched and cached using the neuromaps.datasets module.
By default, data are stored at:

- **macOS/Linux**: `~/neuromaps-data/`

- **Windows**: `C:\Users\<you>\neuromaps-data\`

## References

1. **Markello, R. D., & Mišić, B. (2022).**  
   _Neuromaps: structural and functional maps of the human brain._  
   **Nature Methods, 19**, 1117–1124.  
   [https://doi.org/10.1038/s41592-022-01625-w](https://doi.org/10.1038/s41592-022-01625-w)

2. **Alexander-Bloch, A., Shou, H., Liu, S., Satterthwaite, T. D., Glahn, D. C., Shinohara, R. T., Vandekar, S. N., & Raznahan, A. (2018).**  
   _On testing for spatial correspondence between maps of human brain structure and function._  
   **NeuroImage, 178**, 540–551.  
   [https://doi.org/10.1016/j.neuroimage.2018.05.070](https://doi.org/10.1016/j.neuroimage.2018.05.070)

3. **Hill, J., Inder, T., Neil, J., Dierker, D., Harwell, J., & Van Essen, D. C. (2010).**  
   _Similar patterns of cortical expansion during human development and evolution._  
   **Proceedings of the National Academy of Sciences (PNAS), 107**(29), 13135–13140.  
   [https://doi.org/10.1073/pnas.1001229107](https://doi.org/10.1073/pnas.1001229107)
   _Provdies Evolutionary Expansion source map, and developmental expansion target map._
4. **Vaishnavi, S. N., Vlassenko, A. G., Rundle, M. M., Snyder, A. Z., Mintun, M. A., & Raichle, M. E. (2010).**  
   _Regional aerobic glycolysis in the human brain._  
   **Proceedings of the National Academy of Sciences (PNAS), 107**(41), 17757–17762.  
   [https://doi.org/10.1073/pnas.1010459107](https://doi.org/10.1073/pnas.1010459107)  
   _Covers CBV, CBF, CMRglc, and CMRO₂ metabolic maps used as target datasets._

5. **Glasser, M. F., Sotiropoulos, S. N., Wilson, J. A., Coalson, T. S., Fischl, B., Andersson, J. L., Xu, J., Jbabdi, S., Webster, M., Polimeni, J. R., Van Essen, D. C., & Jenkinson, M. (2013).**  
   _The minimal preprocessing pipelines for the Human Connectome Project._  
   **NeuroImage, 80**, 105–124.  
   [https://doi.org/10.1016/j.neuroimage.2013.04.127](https://doi.org/10.1016/j.neuroimage.2013.04.127)  
   _Reference for HCP S1200-derived myelin and cortical thickness maps._

6. **Margulies, D. S., Ghosh, S. S., Goulas, A., Falkiewicz, M., Huntenburg, J. M., Langs, G., Bezgin, G., Eickhoff, S. B., Castellanos, F. X., Petrides, M., Jefferies, E., & Smallwood, J. (2016).**  
   _Situating the default-mode network along a principal gradient of macroscale cortical organization._  
   **Proceedings of the National Academy of Sciences (PNAS), 113**(44), 12574–12579.  
   [https://doi.org/10.1073/pnas.1608282113](https://doi.org/10.1073/pnas.1608282113)  
   _Provides functional gradient maps 1–3 used in the comparison._

7. **Markello, R. D., Spreng, R. N., Luh, W. M., Anderson, K. M., Ge, T., Holmes, A. J., & Mišić, B. (2021).**  
   _Deconstructing the human connectome: Sensory, transmodal, and association networks are organized along distinct gradients of gene expression._  
   **Nature Neuroscience, 24**, 1255–1265.  
   [https://doi.org/10.1038/s41593-021-00815-8](https://doi.org/10.1038/s41593-021-00815-8)  
   _Source for the Allen Human Brain Atlas gene-expression-derived gradients._

8. **Sydnor, V. J., Larsen, B., Bassett, D. S., Alexander-Bloch, A., Fair, D. A., Liston, C., Mackey, A. P., Milham, M. P., Pines, A., Roalf, D. R., Seidlitz, J., Xu, T., Raznahan, A., Satterthwaite, T. D., & Graham, A. M. (2021).**  
   _Neurodevelopment of the association cortices: Patterns, mechanisms, and implications for psychopathology._  
   **Neuron, 109**(18), 2820–2846.e7.  
   [https://doi.org/10.1016/j.neuron.2021.06.016](https://doi.org/10.1016/j.neuron.2021.06.016)  
   _Associated with the SAaxis map of cortical surface area._

9. **Xu, T., Dong, H.-M., Zhang, S., Zuo, X.-N., & Milham, M. P. (2020).**  
   _Cross-species functional alignment reveals evolutionary hierarchy within the connectome._  
   **NeuroImage, 223**, 117346.  
   [https://doi.org/10.1016/j.neuroimage.2020.117346](https://doi.org/10.1016/j.neuroimage.2020.117346)  
   _Provides Functional Homology (FChomology) and Evolutionary Homology (evoexp) maps._

10. **Reardon, P. K., Seidlitz, J., Vandekar, S., Liu, S., Patel, R., Park, M. T. M., Alexander-Bloch, A., Clasen, L. S., Blumenthal, J. D., Lalonde, F. M., & others. (2018).**  
    _Normative brain size variation and brain shape diversity in humans._  
    **Science, 360**(6394), 1222–1227.  
    [https://doi.org/10.1126/science.aar2578](https://doi.org/10.1126/science.aar2578)  
    _Corresponds to the Allometric Scaling (PNC) and Allometric Scaling (NIH) dataset in CIVET 41k space._
11. **Mueller, S., Wang, D., Fox, M. D., Yeo, B. T. T., Sepulcre, J., Sabuncu, M. R., Shafee, R., Lu, J., & Liu, H. (2013).**  
    _Individual variability in functional connectivity architecture of the human brain._  
    **Neuron, 77**(3), 586–595.  
    [https://doi.org/10.1016/j.neuron.2012.12.028](https://doi.org/10.1016/j.neuron.2012.12.028)  
    _Provides the intersubject functional variability map._
