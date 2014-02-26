\name{ktspshrink}
\alias{ktspshrink}
\title{Shrink a KTSP to given size}
\usage{
ktspshrink(ktsp, dat, k)
}
\arguments{
  \item{ktsp}{a KTSP}

  \item{dat}{training data}

  \item{k}{either an integer specifying the new size of the
  KTSP or a vector specifying wich TSP's should be kept.}
}
\description{
Shrink a KTSP to given size
}
\examples{
data(ktspdata)
ktsp <- ktspcalc(dat,grp,5,"My KTSP")
ktspplotpairs(ktsp,dat)
ktsp.small <- ktspshrink(ktsp,dat,c(1,3,5))
ktspplotpairs(ktsp.small,dat,file="example.pdf")
}

