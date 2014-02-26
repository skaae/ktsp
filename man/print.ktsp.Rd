\name{print.ktsp}
\alias{print.ktsp}
\title{Print KTSP}
\usage{
\method{print}{ktsp}(x, y, ...)
}
\arguments{
  \item{x}{a ktsp}

  \item{y}{unused}

  \item{...}{unused}
}
\description{
Print basic decription of KTSP class
}
\examples{
data(ktspdata)
ktsp <- ktspcalc(dat,grp,5,"Estrogen Receptor")
print(ktsp)
}

