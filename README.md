## Supplementary Materials for Staton and Catalano

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1467684.svg)](https://doi.org/10.5281/zenodo.1467684)

This repository stores the materials for the analysis presented in Staton and Catalano (DOI: 10.1139/cjfas-2018-0176): _Bayesian information updating procedures for Pacific salmon run size indicators: Evaluation in the presence and absence of auxiliary migration timing information_.

### Running this code

To run this code, clone or download this repository to your computer. Open the file `1_Analysis.R` in R, make sure the working directory is set to the location of that file, and run the entire script. The contents of `/Output` will be overwritten when the code finishes in approximately a half hour.

Much of the code is wrapped into functions in the `0_Functions.R` file. For interpretation of this code, we suggest you run through the `1_Analysis.R` code until you run into a user-defined function (`prepare_fit_data()` is the first one). Then sort through the source code for that function until you see how it fits in to the analysis, then move back to the `1_Analysis.R`.
