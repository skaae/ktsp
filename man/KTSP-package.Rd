\docType{package}
\name{KTSP package}
\alias{KTSP package}
\alias{KTSP package-package}
\alias{ktsp}
\title{K-TSP package}
\description{
KTSP package
}
\details{
KTSP is a package for training bundles of Top scoring pair
classifiers. The package include methods for choosing how
the top scoring pairs are ' combined into a single
predction, including agumenting the KTSP with random forest
or neural networks[TODO].

Users are supposed to use \code{ktsp*, print, predict} and
\code{summary} functions \itemize{ \item{\strong{print}
Print name and number of tsps in ktsp object}
\item{\strong{predict} returns vector of predictions}
\item{\strong{summary} prints bacis statistics, e.g.
confusion matrix} \item{\strong{ktspcalc} Train a ktsp with
k tsp's} \item{\strong{ktspbestk} finds optimal number of k
using cross validation} \item{\strong{ktsptrainrf} combines
a ktsp with random forest} \item{\strong{ktsppredmat}
return a matrix of prediction of individual tsp's}
\item{\strong{ktspshrink} remove TSP's from ktsp}
\item{\strong{ktspplotpairs} plot TSP's in ktsp} }

The package build on the switchbox package marchioni et. al
(2013). The original TSP algorihtm was devaloped by Geman
(2004) and extend by Tan (2005)
}
\author{
Soeren Soenderby
}
\references{
Marchionni, Luigi, et al. "A simple and reproducible breast
cancer prognostic test." BMC genomics 14.1 (2013): 336.

D. Geman, C. d'Avignon, D. Naiman and R. Winslow,
"Classifying gene expression profiles from pairwise mRNA
comparisons," Statist. Appl. in Genetics and Molecular
Biology, 3, 2004.

A.C. Tan, D.Q. Naiman, L. Xu, R.L. Winslow, D. Geman,
"Simple decision rules for classifying human cancers from
gene expression profiles," Bioinformatics, 21: 3896-3904,
2005.
}

