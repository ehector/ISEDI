\name{ISEDI-package}
\alias{ISEDI-package}
\alias{ISEDI}
\docType{package}
\title{\packageTitle{ISEDI}}
\description{\packageDescription{ISEDI}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{ISEDI}
  \packageIndices{ISEDI}

  This is a repository for the R package to estimate the parameters of a generalized linear model in a dataset by borrowing information from a prior analysis on another dataset.. The R package's main files are:

- src/KL_funcs.cpp: this file defines the Rcpp functions that compute the likelihoods, KL divergences and their derivatives for the logistic regression model.

-R/KL_funcs.R: this file defines the R function for the ISE estimation of model parameters, as well as functions to compute the mean squared error.
}
\author{
  \packageAuthor{ISEDI}
  
  Maintainer: \packageMaintainer{ISEDI}
}
\references{
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
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
}
\examples{
########### Gaussian example ########### 

## Set parameters for the data generation
n_1 <- 50
n_2 <- 200
p <- 10
beta_1 <- c(1, -1.8, 2.6, 1.4, -3.6, 3.5, 2.4, -3.3, 1.8, -3.4, 2.8)

## Generate the data in the first dataset
set.seed(1234)
X_1 <- cbind(1, matrix(rnorm(n=n_1*p, mean=0, sd=1),ncol=p))
y_1 <- X_1\%*\%beta_1 + rnorm(n=n_1, mean=0, sd=1)
hat_beta_1 <- drop(solve(t(X_1)\%*\%X_1)\%*\%t(X_1)\%*\%y_1)

## Generate the data in the second dataset
set.seed(1357)
beta_2 <- hat_beta_1 + c(0, 0.1, 0, -0.1, 0.1, -0.1, 0, 0, 0, -0.1, 0.1)*2
X_2 <- cbind(1, matrix(rnorm(n=n_2*p, mean=0, sd=1),ncol=p))
y_2 <- X_2\%*\%beta_2 + rnorm(n=n_2, mean=0, sd=1)

## Obtain the ISE estimates
lambda_seq <- seq(0,10,length.out=1000)
estimates <- LDW_meta(X_1, X_2, y_1, y_2, lambda_seq, family="gaussian")

## Optimal estimates are returned in:
estimates$estimates_opt

########### Bernoulli example ########### 

## Set the parameters for the data generation
n_1 <- 500
n_2 <- 500
p <- 4
beta_1 <- c(1, -1.8, -1.2, 1.6, 0.2)

## Generate the data in the first dataset
set.seed(1234)
X_1 <- cbind(1, matrix(rnorm(n=n_1*p, mean=0, sd=0.75),ncol=p))
eta_1 <- exp(X_1\%*\%beta_1)/(1+exp(X_1\%*\%beta_1))
y_1 <- rbinom(n=n_1, size=1, prob=eta_1)
hat_beta_1 <- c(coef(glm(y_1 ~ 0 + X_1, family="binomial")))

## Generate the data in the second dataset
set.seed(1357)
beta_2 <- hat_beta_1 + rep(0,5)
X_2 <- cbind(1, matrix(rnorm(n=n_2*p, mean=0, sd=1),ncol=p))
eta_2 <- exp(X_2\%*\%beta_2)/(1+exp(X_2\%*\%beta_2))
y_2 <- rbinom(n=n_2, size=1, prob=eta_2)

## Obtain the ISE estimates
lambda_seq <- seq(0,7,length.out=100)
estimates <- LDW_meta(X_1, X_2, y_1, y_2, lambda_seq, family="binomial")

## Optimal estimates are returned in:
estimates$estimates_opt
}
