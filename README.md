# ISEDI
Information Sharing for Efficient Data Integration

This is a repository for the R package to estimate the parameters of a generalized linear model in a dataset by borrowing information from a prior analysis on another dataset. The R package's main files are:
- src/KL_funcs.cpp: this file defines the Rcpp functions that compute the likelihoods, KL divergences and their derivatives for the logistic regression model.
- R/KL_funcs.R: this file defines the R function for the ISE estimation of model parameters, as well as functions to compute the mean squared error.

The ISEDI man file contains an example for running the regression models from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The ISEDI R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL ISEDI_1.0-1.tar.gz
- from the downloaded and renamed ISEDI folder as R CMD build ISEDI and R CMD INSTALL ISEDI_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the ISEDI R package, please consider citing the relevant manuscript: E.C. Hector and R. Martin (2022+). Turning the information-sharing dial: efficient inference from different data sources. arXiv, arXiv:2207.08886.

# References

Efron, B. and Morris, C. (1977). Stein’s paradox in statistics. Scientific American, 236(5):119–127.

Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased estimation for nonorthogonal problems. Technometrics, 12(1):55–67.
