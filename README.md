# Code to reproduce the BIOCARD analysis in "Bayesian shape-constrained regression for quantifying Alzheimer's disease progression with biomarkers" by Li et al. 2024

This repository contains instructions and scripts used to reproduce the BIOCARD analysis in *Bayesian shape-constrained regression for quantifying Alzheimer's disease progression with biomarkers* by Li et al. 2024.

## Setting up Dependencies

The scripts uses multiple existing packages in R, including **tidyverse, splines2, TruncatedNormal, mvtnorm, matrixStats, foreach, progress, doParallel, doRNG, VGAM, abind** and **hdtg**. The first nine packages can be installed directly through CRAN by using the following code in R console:

```         
install.packages(c(
'tidyverse',
'splines2',
'TruncatedNormal',
'mvtnorm',
'matrixStats',
'foreach',
'progress',
'doParallel',
'doRNG',
'VGAM',
'abind'
))
```

The **hdtg** package can be installed through GitHub using the following code in R console:

```         
remotes::install_github("https://github.com/suchard-group/hdtg", build = FALSE)
```

## Reproducing the Results in Section 5: Real Data Application

This repository contains three R script files used to analyze BIOCARD study data and produce graph and table results. To reproduce the results in Section 5, run the following code with working directory set to the one containing the three R files, and with the data files *biocard.rda, AD_Diagnose_Age.rda, DECAGE.rda* placed in the same directory of these scripts:

```         
Rscript main_S_realdata.R
Rscript manuscript_plot.R
```

After running these scripts, there should be one text file and two PDF files. *table.txt* contains **Table 2** in Section 5.4, which depicts the point estimate and credible intervals of effects of age and other covariates. Figures 2 and 3 are contained in *SplineStd.pdf* and *Goodness_of_fit_Grid3.pdf*, which depicts Bayes estimates of standardized biomarker progression curves, and spaghetti plots of observed and fitted outcomes for select biomarkers.

There should also be a workspace file *biocard_result_group20nonzeros.RData* after running the scripts. It contains all the MCMC samples for reference.
