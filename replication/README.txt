Data and codes for "Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes" by Chernozhukov, Fernandez-Val, Wuethrich and Melly

This folder contains:

1) The R package discreteQ
This package implements the procedure suggested in the paper. The package can be installed
from the source file "discreteQ_0.3.0.tar". The latest updated version of the package can
be installed by typing
> devtools::install_github("bmelly/discreteQ", build_vignettes = TRUE)
in R. After installing and loading the package, we recommend reading the help file for the command discreteQ and the vignette:
> library(discreteQ)
> ?discreteQ
> vignette("discreteQ", package="discreteQ")

2) The datasets used in the applications
The folder "data" contains the datasets "data_oregon.dta" (used in the first application
in Section 4.3) and "cpp_selected.dta" (used in the second application in Section 4.4).

3) The R codes that replicates all the figures and tables of the paper
The folder "codes" contains the R codes that replicate all the figures and tables in the
paper. The path to the folder containing the data and codes must be modified in the
file "_master.R". This file will then call the other documents, which do not need to be
modified. Each document produces the figure or table indicated in its name. Note that some
results take a very long time to be computed (up to 10 days for the second application)
because the number of observations and covariates is large. The approximate computing
time are given in _master.R for each figure.
