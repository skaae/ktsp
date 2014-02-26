\name{predict.ktsp}
\alias{predict.ktsp}
\title{Calculate prediction from KTSP object}
\usage{
\method{predict}{ktsp}(object, data, combineFunc = majorityPredict, ...)
}
\arguments{
  \item{object}{A KTSP object}

  \item{data}{a matrix with rows as genes and samples as
  columns}

  \item{combineFunc}{A function for combining the k-tsp's
  into a single prediction. Default is majority predict.
  See \code{\link{majorityPredict}} for example}

  \item{...}{not used}
}
\value{
vector of predictions
}
\description{
Takes a KTSP and data as input and returns the predictions
as a vector. Optionally takes a function for combining
output from the k-tsp's into a single prediction for each
sample. If no combine function is given majority prediction
is used.
}
\examples{
data(ktspdata)
ktsp <- ktspcalc(dat,grp,5)
pred <- predict(ktsp,dat)
}

