\name{cvevalrf}
\alias{cvevalrf}
\title{Evaluate random forest performance on cv splits}
\usage{
cvevalrf(dat, grp, lst.ktsp, mtry, nodesize, features, evalfun)
}
\arguments{
  \item{lst.ktsp}{list of ktsp, format as output of
  \code{\link{cvtrain}}}

  \item{dat}{gene expression matrix. Genes as rows and
  samples as columns}

  \item{grp}{binary encoded labels}

  \item{evalfun}{evaluation function see
  \code{\link{CalcAUC}}}

  \item{mtry}{mtry hyperparam vals}

  \item{nodesize}{nodesize hyper param vals}

  \item{features}{hyperparam vals}
}
\value{
df with the performance for TSP's up to k. Columns: #'
\item{Train}{Training performance} \item{Test}{Test
performance} \item{k}{Number of included tsp's}
}
\description{
Evaluate random forest performance on cv splits
}

