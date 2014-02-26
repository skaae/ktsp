\name{cvsplits}
\alias{cvsplits}
\title{Create Crossvalidation split indices}
\usage{
cvsplits(Nobs, K = 5, seed = 1234)
}
\arguments{
  \item{Nobs}{Number of observations in dataset}

  \item{K}{fold cross validation}

  \item{seed}{seed for cross validation splits}
}
\value{
a list of length ??code(k) with train and test indices
}
\description{
Create Crossvalidation split indices
}

