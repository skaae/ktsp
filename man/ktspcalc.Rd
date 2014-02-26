\name{ktspcalc}
\alias{ktspcalc}
\title{train KTSP classifier.}
\usage{
ktspcalc(data, grp, n, name = NULL)
}
\arguments{
  \item{data}{gene matrix. Rows are genes, columns are
  samples}

  \item{grp}{a vector of training labels encoded as 1 and
  0}

  \item{n}{The number of TSP's to include.}

  \item{name}{optional name of KTSP. For plots etc}
}
\value{
Retuns a KTSP object with the following fields
\item{TSPs}{Matrix of indices of selected pairs. k by 2}
\item{k}{Number of TSP's} \item{geneNames}{Matrix of names
of selected pairs. k by 2} \item{name}{Optional name of
ktsp used for plotting etc} \item{train.grp}{labels of
training data} \item{train.pred}{prediction on training
data}
}
\description{
Train a KTSP classifier with n TSP's.
}
\examples{
data(ktspdata)
ktsp <- ktspcalc(dat,grp,5,"My KTSP")
summary(ktsp)
}
\keyword{ktsp}
\keyword{train}

