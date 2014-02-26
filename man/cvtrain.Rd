\name{cvtrain}
\alias{cvtrain}
\title{Train a ktsp on each cv split}
\usage{
cvtrain(splits, dat, grp, k, name)
}
\arguments{
  \item{splits}{cv splits see \code{\link{cvsplits}}}

  \item{dat}{gene expression matrix. Genes as rows and
  samples as columns}

  \item{grp}{binary encoded labels}

  \item{k}{number of tsp's}

  \item{name}{TSP name}
}
\value{
list of length "number of corssval folds" with the
following fields \item{ktsp}{a ktsp object}
\item{idx.train}{indices of cv train data}
\item{idx.test}{indices of cv test data}
}
\description{
This function is used to evaluate CV performance e.g. in
\code{\link{ktspbestk}} and \code{\link{ktsptrainrf}}
}

