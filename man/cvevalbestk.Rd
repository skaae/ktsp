\name{cvevalbestk}
\alias{cvevalbestk}
\title{Evaluate performance seq search CV}
\usage{
cvevalbestk(dat, grp, lst.ktsp, evalfun, reg, k.max)
}
\arguments{
  \item{lst.ktsp}{list of ktsp, format as output of
  \code{\link{cvtrain}}}

  \item{dat}{gene expression matrix. Genes as rows and
  samples as columns}

  \item{grp}{binary encoded labels}

  \item{evalfun}{evaluation function see
  \code{\link{CalcAUC}} for example}
}
\value{
df with the performance for TSP's up to k. Columns:
\item{Train}{Training performance} \item{Test}{Test
performance} \item{k}{Number of included tsp's}
}
\description{
Evaluate performance seq search CV
}

