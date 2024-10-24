\name{mixSVG2}
\alias{mixSVG2}
\title{mixSVG: An Omnibus Test of Mixed Effects for Detecting Spatially Variable Genes with Spatial Transcriptomics Data}
\usage{
mixSVG2(count,
       coord,
       X = NULL,
       libsize_inc = TRUE,
       libsize = NULL,
       vtest_zero_prop = 0.995,
       ncore = 10,
       n_perm = 1000,
       sig = 0.05)
}
\description{
The mixSVG function implements an omnibus test to detect spatially variable genes with spatial transcriptomics data, which tests for fixed effect (mean-level) and random effect (variance-level) of spatial heterogeneity of gene expression simultaneously.
}

\arguments{
\item{count}{A q by n numeric matrix of gene expression counts for q genes and n spots. }
\item{coord}{An n by 2 numeric matrix of two-dimensional spatial coordinates for n spots.}
\item{X}{A n by p numeric matrix of p covariates for n spots.}
\item{libsize_inc}{Whether to account for the library size. \code{libsize_inc=TRUE} by default.}
\item{libsize}{A numeric vector of the library size for n spots. If \code{libsize_inc=TRUE}, then \code{libsize} will be the total expression counts of all genes on each spot by default. If \code{libsize_inc=FALSE}, then \code{libsize=1} by default.}
\item{vtest_zero_prop}{A numeric value between 0 and 1. The mixed effects (fixed effect and random effect) will be tested when the proportion of zero expression counts across all spots is less than \code{vtest_zero_prop}. Otherwise, only the fixed effect will be tested.}
\item{ncore}{Number of cores used for parallel computation}
\item{n_perm}{Number of permutation replicates for testing the random effect.}
\item{sig}{Significance level for detecting spatially variable genes based on the adjusted P-values of Benjamini-Hochberg method.}
}

\value{
\item{results}{A list of results for all genes, including:
\itemize{
\item \code{model0}: The estimation results under the null model, which contain the estimated coefficients of covariates (\code{beta}), the estimated variance of residual (\code{tau}), an indicator of the estimation convergence (\code{converge}), the number of iterations (\code{iter});
\item \code{pval}: The combined P-value of mixSVG based on 13 kinds of transformations of spatial coordinates accounting for different spatial patterns;
\item  \code{pval_pat}: A 13 by 3 matrix of P-values of mixSVG, where 13 rows represent 13 transformations of spatial coordinates, and three columns represent the test for mixed effect, fixed effect and random effect, respectively.
}
}

\item{pval_all}{A matrix of P-values for all genes. The first column contains the original P-values of mixSVG, and the second column contains the P-values adjusted by the Benjamini-Hochberg method.}

\item{pval_sig}{A matrix of P-values for the detected spatially variable genes (whose adjusted P-values > \code{sig}). The first column contains the original P-values of mixSVG, and the second column contains the P-values adjusted by the Benjamini-Hochberg method.}
}









