# Code to reproduce the BIOCARD analysis in "Bayesian shape-constrained regression for quantifying Alzheimer's disease progression with biomarkers" by Li et al. 2024

This repository contains instructions and scripts used to reproduce the BIOCARD analysis in *Bayesian shape-constrained regression for quantifying Alzheimer's disease progression with biomarkers* by Li et al. 2024.

## Synthetic Data

This repository contains the three *rda* files needed for the analysis, but the data within are synthetic and resembles the actual preprocessed data in format only. These files are solely for **demonstration purposes** and should **NOT** be used or regarded as actual BIOCARD data. The following table contains the codebook for these files.

| Variable              | Description | Note |
| :---------: | :----------------------: | :-----------------: |
|*biocard.rda*|||
|MMSCORE|Mini Mental State Examination|NA allowed|
|logmem|Logical Memory Test|NA allowed|
|DSST|Digit Symbol Substitution Test|NA allowed|
|biec.thik|Entorhinal cortex thickness|NA allowed|
|Hippo_dadjust|Hippocampus volume|NA allowed|
|Ent_dadjust|Entorhinal cortex volume|NA allowed|
|MTL1|MRI standardized residual composites|NA allowed|
|SPARE_AD|Spatial Pattern of Abnormality for Recognition of Early Alzheimer’s Disease score|NA allowed|
|ttau|total tau|NA allowed|
|ptau181|phosphorylated tau|NA allowed|
|AB42AB40| $A\beta_{1-42}$ to $A\beta_{1-40}$ ratio|NA allowed|
|apoe|APOE4 carrier status|Carrier = 1, Non-carrier = 0|
|SEX|biological sex|Female = 1, Male = 0|
|education|(standardized) years of education||
|Study.ID|Individual ID||
|ageori|Age (in years)||
|*DECAGE.rda*|||
|decs|Symptom onset age||
|*AD_Diagnose_Age.rda*|||
|AD_Diagnose_Age|Dementia onset age||


## Setting up Dependencies

The scripts uses multiple existing packages in R, including **tidyverse, splines2, TruncatedNormal, mvtnorm, matrixStats, foreach, progress, doParallel, doRNG, VGAM, abind** and **hdtg**. The first nine packages can be installed directly through CRAN by using the following code in R console:

```         
install.packages(c(
'latex2exp'
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

