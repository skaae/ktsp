\name{ktsppredmat}
\alias{ktsppredmat}
\title{Calculate prediction from ktsp}
\usage{
ktsppredmat(data, ktsp)
}
\arguments{
  \item{data}{expresssion matrix. Genes are rows and
  samples are columns}

  \item{ktsp}{ktsp}
}
\value{
matrix of predictions. size samples by number_of_tsp
}
\description{
Returns a matrix of size samples by number_of_tsp Each row
is the TSP prediction of for a sample
}
\examples{
#Example 1
data(ktspdata)
ktsp <- ktspcalc(dat,grp,5,name="Test TSP")
predmat <- ktsppredmat(dat,ktsp)
pred <- apply(predmat,2,majorityPredict)
AUC <- CalcAUC(pred,ktsp$train.grp)
stopifnot(AUC==CalcAUC(ktsp$train.pred,ktsp$train.grp))

#Example 2 - load data and precalculated tsp
data(her2)
data(ktsp)
predmat <- ktsppredmat(expr.train.t.knn,ktsp)
}

