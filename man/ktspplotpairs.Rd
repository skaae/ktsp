\name{ktspplotpairs}
\alias{ktspplotpairs}
\title{Plot KTSP pairs}
\usage{
ktspplotpairs(ktsp, dat, file = NULL)
}
\arguments{
  \item{ktsp}{a KTSP}

  \item{dat}{training data}

  \item{file}{optionally set save path.}
}
\description{
Plots all pairs in a TSP. Use \code{\link{ktspshrink}} if
you do only want to plot some of the pairs
}
\examples{
data(ktspdata)
ktsp <- ktspcalc(dat,grp,5,"My KTSP")
ktspplotpairs(ktsp,dat)
ktspplotpairs(ktsp,dat,file="example.pdf")
}

