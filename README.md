# 2DSFS-scan
Python scripts for inference of genomic divergence between populations using the joint site frequency spectrum (SFS) and composite likelihoods. The approach aims to detect aberrant 2D SFS relative to the background genome. The pipeline uses VCF and population map files as input to extract SNP data, computes the 1D and 2D SFS, and calculate likelihood statistics in genomic windows of a fixed size.

## Code

To compile our code for calculating T2D, run `python setup.py` from the top-level directory. The main source code is contained in `src/twoDSFS/likelihoodCalcs.py`. Simulation scripts are found in `slim_sims/stabSelGrad_model/simulations_LHU/StabSel_dadiparams_LHU.slim` and `slim_sims/stabSelGrad_model/simulations_LHU/StabSel_dadiparams_LHU.noDel.slim`.

## Data

VCF files for ECB empirical data are too large for GitHub and will be uploaded to a separate repository upon publication.

## Figures

### Figure 1
To generate Figure 1, run the R script `Figures/`

### Figure 2
To generate Figure 2, run the R script `Figures/Fig1SFS.R` 

### Figures 3 and 4
To generate Figures 3 and 4, run the R script `Figures/power.R`

## Figure 5
To generate Figure 5, run the R script `Figures/`

## Figure 6
To generate Figure 6, run the R script `Figures/fst_T2D_emp.R`


