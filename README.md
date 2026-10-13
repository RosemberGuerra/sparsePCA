# A Guide for Sparse PCA: Model Comparison and Applications

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![DOI](https://img.shields.io/badge/DOI-10.1007%2Fs11336--021--09773--2-blue)](https://doi.org/10.1007/s11336-021-09773-2)

> **Replication code and comparative benchmarks across Sparse PCA formulations, implementations, and empirical applications.**

This repository contains the complete MATLAB and R code supporting the methodological review and empirical benchmark published in:

> Guerra-Urzola, R., Van Deun, K., Vera, J. C., & Sijtsma, K. (2021). **A Guide for Sparse PCA: Model Comparison and Applications**. *Psychometrika*, 86(4), 893–919. [https://doi.org/10.1007/s11336-021-09773-2](https://doi.org/10.1007/s11336-021-09773-2)

---

## Repository Structure

```text
├── Applications/
│   ├── big5/                     # Big Five personality inventory benchmark
│   └── gene_expression/          # Microarray gene expression_application
├── data_generation/               # Synthetic matrix generation (loading & score structures)
├── Numerical_Experiments/
│   ├── R/                         # R implementations (Varimax, Simplimax, SPCA)
│   ├── matlab/                    # MATLAB implementations (GPower, pathSPCA, sPCArSVD)
│   └── Plots/                     # Visualizations and trade-off comparison scripts
└── README.md
```
## External Dependencies
### R Dependencies
Requires R (>= 4.0.0). Install CRAN dependencies:
```R
install.packages(c("elasticnet", "GPArotation", "psych", "ggplot2", "dplyr"))
```
### MATLAB Dependencies
Most functions are self-contained. Two external solvers require downloading original author utilities:
- GP.m / GPower: Requires GPower.m from Michel Journée's website. Place GPower.m in the working directory or MATLAB path.
- pathSPCA.m: Requires maxeig.m and FullPathCov.m from Alexandre d'Aspremont's website. Place both files in the working directory or MATLAB path.
## Reproduction Workflows
### 1. Data Generation
To generate synthetic benchmark datasets with sparse loading patterns ($P$ and $W$):
- In R: Run scripts inside Data Generation/ (e.g., PandWsparse.R, Psparse.R, Wsparse.R).
- In MATLAB: Run corresponding .m scripts (e.g., PandWsparse.m).
### 2. Numerical Simulation Studies
To evaluate recovery error, variance explained, and sparsity trade-offs:
- R simulations:
```R
source("Numerical Experiments/R/SPCA_simulation.R")
source("Numerical Experiments/R/VARIMAX_simulation.R")
source("Numerical Experiments/R/SIMPLIMAX_simulation.R")
```
- MATLAB simulations:
```Matlab
% Run from MATLAB command line within Numerical Experiments/matlab/
run('GPower_simulation.m')
run('pathSPCA_simulation.m')
run('sPCArSVD_simulation.m')
```
- Plotting results:
```R
source("Numerical Experiments/Plots/JoiningDataSets.R")
source("Numerical Experiments/Plots/PlotSimulations.R")
```
### 3. Empirical Applications
#### Big Five Personality Data
1. Download and preprocess data:
```R
source("Applications/big5/Download&Processing.R")
```
2. Run comparative models (R and MATLAB scripts in Applications/big5/).
3. Generate comparison figures:
```R
source("Applications/big5/Plots.R")
```
#### Gene Expression Data
1. Download and process microarray dataset:
```R
source("Applications/Gene Expression/Download&Processing.R")
```
2. Run Sparse PCA formulations (Application_GPower.m, Application_spcarsvd.m, ISvsPSvsPEV.R).

## Citation
If you use this benchmark suite or implementations in your research, please cite the original article:
```Code snippet
@article{guerra_urzola_guide_2021,
  author    = {Guerra-Urzola, Rosember and Van Deun, Katrijn and Vera, Juan C. and Sijtsma, Klaas},
  title     = {A Guide for Sparse PCA: Model Comparison and Applications},
  journal   = {Psychometrika},
  volume    = {86},
  number    = {4},
  pages     = {893--919},
  year      = {2021},
  doi       = {10.1007/s11336-021-09773-2}
}
```
## License
This repository is distributed under the MIT License - see the LICENSE file for details.