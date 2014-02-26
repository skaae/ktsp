\name{evalfs}
\alias{evalfs}
\title{Select bestk using forward/backward feature selection}
\usage{
evalfs(attrib, predmat, grp, evalfun)
}
\arguments{
  \item{lst.ktsp}{list of ktsp, format as output of
  \code{\link{cvtrain}}}

  \item{dat}{gene expression matrix. Genes as rows and
  samples as columns}

  \item{grp}{binary encoded labels}

  \item{evalfun}{evaluation function see
  \code{\link{CalcAUC}} for example}

  \item{type}{"forward" or "backward" search}

  \item{predmats}{list of precalculated prediction
  matrices}

  \item{dat}{gene expression df Genes as rows and samples
  as columns}

  \item{grp}{binary encoded labels}

  \item{cross}{number of cross validation splits}

  \item{attrib}{features to try}

  \item{evalfun}{evaluation function see
  \code{\link{CalcAUC}} for example}
}
\value{
list with the following fields
\item{combineFunc}{CombineFunc to be used with a ktsp
trained with k tsp's} \item{k}{train ktsp with the many
features} \item{features}{selected features} evalute
performance, used in forward/backward search

list with the following fields
\item{combineFunc}{CombineFunc to be used with a ktsp
trained with k tsp's} \item{k}{train ktsp with the many
features} \item{features}{selected features}
}
\description{
Select bestk using forward/backward feature selection
}

