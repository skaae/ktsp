\name{ktsptrainnn}
\alias{ktsptrainnn}
\title{Create Neural net combine function}
\usage{
ktsptrainnn(dat, grp, k.max, cross, evalfun = CalcAUC, splits = NULL,
  seed = 1234, size = NULL, decay = NULL, features = NULL,
  precalc = NULL)
}
\arguments{
  \item{dat}{gene expression matrix. Genes as rows and
  samples as columns}

  \item{grp}{binary encoded labels}

  \item{k.max}{largest TSP to evaluate}

  \item{cross}{number of crosss validation folds}

  \item{evalfun}{function for eval performance on CV
  splits. Defaults to AUC}

  \item{splits}{optionally supply the cv split indices}

  \item{seed}{seed for cross validation splits}

  \item{decay}{neuralnet weight penalty, vector for several
  values}

  \item{size}{hiddenlayer sizes to try}

  \item{features}{vector of number of features to try}
}
\value{
A list with the following fields \item{hyperparams}{A
matrix with average test/trainset performance using
different hyperparameter values} \item{bestparams}{best
hyper parameters } \item{ktsp}{A ktsp trained on the
training data using best k tsp's} \item{combineFunc}{The
combine function used/created}
}
\description{
Using CV this function trains a Neural Net ontop of a KTSP.
User specify hyper parameter space to be searched.
}
\examples{
\dontrun{
data(ktspdata)
cross = 5
size =
decay = c(5e-02, 5e-03, 5e-04, 5e-05)
features = c(5, 10, 25, 50)
splits <- cvsplits(ncol(dat), cross)
res.nn <- ktsptrainnn(dat,grp,100,cross,splits=splits,
                   size =c(20,30,40,50),
                   decay=c(5e-02, 5e-03, 5e-04, 5e-05),
                   features=c(5, 10, 25, 50))
pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)

require(datathesis)
data(her2)
cross = 5
res.her2.nn <- ktsptrainnn(dat,grp,200,cross,splits=splits,
                   size =c(40,50,75,100),
                   decay=c(5,5e-01,5e-02, 5e-03, 5e-04, 5e-05),
                   features=c(5,7,10, 25, 50))
pred.train <- predict(res.her2.nn$ktsp,expr.train.t.knn,
                combineFunc=res.her2.nn$combineFunc)
auc.train <- CalcAUC(pred,clin.train$her2)
summary(besttsp)
}
}
\seealso{
Other hyperparams: \code{\link{ktspbestk}};
\code{\link{ktsptrainrf}}
}

