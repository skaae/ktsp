\name{ktspbestk}
\alias{ktspbestk}
\title{Determine best K using cross validation}
\usage{
ktspbestk(dat, grp, k.max, cross, type, evalfun = CalcAUC, name = NULL,
  splits = NULL, file = NULL, seed = 1234, precalc = NULL,
  reg = noReg)
}
\arguments{
  \item{dat}{gene expression matrix. Genes as rows and
  samples as columns}

  \item{grp}{binary encoded labels, vector}

  \item{k.max}{largest KTSP to evaluate, may be omitted if
  precalc is supplied}

  \item{cross}{number of crosss validation folds,may be
  omitted if precalc is supplied}

  \item{file}{optional file name for plot. Ignored if fs
  =TRUE}

  \item{seed}{seed for cross validation splits}

  \item{evalfun}{function for evaluating performance on CV
  splits. Defaults to AUC}

  \item{name}{Optional name for plot}

  \item{splits}{optionally supply the cv split indices, if
  supplied cross and seed is ignored}

  \item{type}{Feature selection type. If "seq" the function
  tries TSP's sequentially, i.e starting with a ktsp with
  TSP-1, then a KTSP with TSP-1+2, 1+2+3 etc. If "forward"
  greedy forward search, if "backward", greedy backward
  search}

  \item{precalc}{optionally supply a precalculated list of
  CV tsps, output from \code{\link{cvtrain}}, to reduce
  computational time. If supplied k.max, cross and seed are
  ignored.}

  \item{reg}{regularization function. see
  \code{\link{Reg}}}
}
\value{
A list with the following fields \item{features}{either
bestk value (fs = "seq") or the selcted features (fs =
"backward/forward")} \item{ktsp}{A ktsp trained on the
training data using best k tsp's} \item{combineFunc}{The
combine function used/created}
}
\description{
Calculates the optimal number of TSP's using cross
validation
}
\examples{
\dontrun{
# Example 1 - selecting top TSP's (no FS)
data(ktspdata)
res <- ktspbestk(dat,grp,20,3,file="example.pdf",name="Test TSP",
                type="seq")
pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
summary(res$ktsp)

# forward search
res <- ktspbestk(dat,grp,20,3,file="example.pdf",name="Test TSP",
                type="forward")
pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)

# backward search
res <- ktspbestk(dat,grp,20,3,file="example.pdf",name="Test TSP",
                type="backward")
pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)

# Example - her2
require(datathesis)
data(her2) # load expression/clinical data and precalculated ktsp list
resFSseq<- ktspbestk(expr.train.t.knn,clin.train$her2,type="seq",
                       file="example.pdf",name="HER2 FS:sequential",
                       precalc=lst.ktsp)
resFSf <- ktspbestk(expr.train.t.knn,clin.train$her2,type="forward",
                       file="example.pdf",name="HER2 FS:Forward",
                       precalc=lst.ktsp)
resFSb <- ktspbestk(expr.train.t.knn,clin.train$her2,type="backward",
                       file="example.pdf",name="HER2 FS:Backward",
                       precalc=lst.ktsp)

predFSseq  <- predict(resFSseq$ktsp,expr.test.t.knn,
                        combineFunc=resFSseq$combineFunc)
predFSf  <- predict(resFSf$ktsp,expr.test.t.knn,
                         combineFunc=resFSf$combineFunc)
predFSf  <- predict(resFSb$ktsp,expr.test.t.knn,
                         combineFunc=resFSb$combineFunc)
aucFSseq <- CalcAUC(predFSseq,clin.test$her2)
aucFSf  <- CalcAUC(predFSf,clin.test$her2)
aucFSb  <- CalcAUC(predFSb,clin.test$her2)
}
}
\seealso{
Other hyperparams: \code{\link{ktsptrainnn}};
\code{\link{ktsptrainrf}}
}

