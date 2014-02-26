\name{ktsptrainrf}
\alias{ktsptrainrf}
\title{Create RandomForest combine function}
\usage{
ktsptrainrf(dat, grp, k.max, cross, evalfun = CalcAUC, splits = NULL,
  seed = 1234, mtry = NULL, nodesize = NULL, features = NULL,
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

  \item{mtry}{randomforest hyperparam
  \code{link{randomForest}}. Specify values to be tried as
  vector}

  \item{nodesize}{randomfforest hyperparam
  \code{link{randomForest}}. Specify values to be tried as
  vector}

  \item{features}{vector of number of features to try see
  \code{\link{CalcAUC}}}
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
Using CV this function trains a RandomForest ontop of a
KTSP. User specify hyper parameter space to be searched.
}
\examples{
data(ktspdata)
cross = 5
splits <- cvsplits(ncol(dat), cross)
res <- ktsptrainrf(dat,grp,100,cross,splits=splits, mtry=c(1,2,3),
                   nodesize=c(1,2,3,4,5),features=c(10,50,100))
pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
summary(besttsp)

data(her2)
cross = 5
k.max = 200
res <- ktsptrainrf(expr.train.t.knn,clin$her2,k.max,cross,splits=splits,
                   mtry=c(1,2,3),nodesize=c(1,2,3,4,5),
                   features=c(5,10,30,50,100,200))
pred <- predict(res$ktsp,expr.train.t.knn,combineFunc=res$combineFunc)
auc <- CalcAUC(res$ktsp,clin.test$her2)
}
\seealso{
Other hyperparams: \code{\link{ktspbestk}};
\code{\link{ktsptrainnn}}
}

