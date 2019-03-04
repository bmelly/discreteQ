
<!-- README.md is generated from README.Rmd. Please edit that file -->
The discreteQ repository
========================

This repository contains two different parts:

1.  In the folder repository you can find the datasets and the R codes that generated all the results in the paper in the paper "Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes" written by Victor Chernozhukov, Ivan Fernandez-Val, Blaise Melly and Kaspar Wuethrich. This paper is available at <https://arxiv.org/abs/1608.05142>. The first empirical example uses data from the Oregon health experiment. The data are available at <http://www.nber.org/oregon/4.data.html>. The second empirical example is a decomposition of the black-white testscore gap. The data are available at dx.doi.org/10.1257/aer.103.2.981. The codes that produce the simulations results in the supplementary appendix are also provided.

2.  The R package that contain generic functions that allows researchers to easily apply the methods suggested in the paper "Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes" written by V. Chernozhukov, I. Fernandez-Val, B. Melly, and K. Wuethrich. This paper is available at <https://arxiv.org/abs/1608.05142>. The goal of discreteQ is to perform inference on quantile functions, quantile treatment effect functions and decompositions of differences between quantile functions for possibly discrete outcomes.

Installation of the discreteQ package
=====================================

You can install the released version of discreteQ from [CRAN](https://CRAN.R-project.org) with:

``` r
# You can get the discreteQ package from GitHub:
install.packages("devtools")
devtools::install_github("bmelly/discreteQ", build_vignettes = TRUE)
```

After installing and loading the package, we recommend reading the help file for the command discreteQ and the vignette: &gt; library(discreteQ) &gt; ?discreteQ &gt; vignette("discreteQ", package="discreteQ")

``` r
library(discreteQ)
vignette("discreteQ", package="discreteQ")
?discreteQ
```

Concluding remarks
==================

The specificty of discreteQ is that it can deal with discrete and mixed continuous-discrete outcomes. For continuous outcomes we recommend to use the package Counterfactual that can be installed from CRAN.

For the precise definition of the estimators and their statistical properties.You are invited to read the paper "Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes" written by Victor Chernozhukov, Ivan Fernandez-Val, Blaise Melly, and Kaspar Wuethrich. This paper is available at <https://arxiv.org/abs/1608.05142>.
