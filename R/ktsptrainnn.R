#' Create Neural net combine function
#' 
#' Using CV this function trains a Neural Net ontop of a KTSP. User specify 
#' hyper parameter space to be searched.
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param k.max largest TSP to evaluate
#' @param cross number of crosss validation folds 
#' @param evalfun function for eval performance on CV splits. Defaults to AUC 
#' @param splits optionally supply the cv split indices
#' @param seed seed for cross validation splits
#' @param decay neuralnet weight penalty, vector for several values
#' @param size hiddenlayer sizes to try
#' @param features vector of number of features to try
#' @export
#' @import nnet
#' @return A list with the following fields
#'      \item{hyperparams}{A matrix with average test/trainset performance using
#'      different hyperparameter values}
#'      \item{bestparams}{best hyper parameters }
#'      \item{ktsp}{A ktsp trained on the training data using best k tsp's}
#'      \item{combineFunc}{The combine function used/created}
#' @family hyperparams
#' @examples
#' \dontrun{
#' data(ktspdata)
#' cross = 5
#' size = 
#' decay = c(5e-02, 5e-03, 5e-04, 5e-05)
#' features = c(5, 10, 25, 50) 
#' splits <- cvsplits(ncol(dat), cross)
#' res.nn <- ktsptrainnn(dat,grp,100,cross,splits=splits, 
#'                    size =c(20,30,40,50),
#'                    decay=c(5e-02, 5e-03, 5e-04, 5e-05),
#'                    features=c(5, 10, 25, 50))
#' pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
#' 
#' require(datathesis)
#' data(her2)
#' cross = 5
#' res.her2.nn <- ktsptrainnn(dat,grp,200,cross,splits=splits, 
#'                    size =c(40,50,75,100),
#'                    decay=c(5,5e-01,5e-02, 5e-03, 5e-04, 5e-05),
#'                    features=c(5,7,10, 25, 50))
#' pred.train <- predict(res.her2.nn$ktsp,expr.train.t.knn,
#'                 combineFunc=res.her2.nn$combineFunc)
#' auc.train <- CalcAUC(pred,clin.train$her2)
#' summary(besttsp)
#' }
ktsptrainnn <- function(dat, grp, k.max, cross, evalfun = CalcAUC,
                        splits = NULL, seed = 1234,size=NULL,decay=NULL,
                        features=NULL,precalc=NULL){
    stopifnot(all(!c(is.null(decay),is.null(size),is.null(features))==T))
    stopifnot(nrow(grp) == nrow(dat))
    if(is.null(precalc)){
        if(is.null(splits)){
            splits   <- cvsplits(ncol(dat),cross,seed)
        }
        lst.ktsp <- cvtrain(splits,dat,grp,k.max,name)
    }else{
        lst.ktsp <- precalc
        cross = length(lst.ktsp)
        k.max <- lst.ktsp[[1]]$ktsp$k
    }
    
    pf       <- cvevalnn(dat,grp,lst.ktsp,evalfun,size=size,features=features,
                         decay=decay)
    besthyper <- pf[which.max(pf$Test), ,drop=F]
    nn.models <- trainbestnn(besthyper,lst.ktsp,dat,grp)
    
    ktsp      <- ktspcalc(dat,grp,k.max,"test")
    ping(sound=2)
    
    predmat   <- t(ktsppredmat(dat,ktsp))
    
    model <- nnet(
        predmat[ ,seq_len(besthyper$features),drop=F], 
        grp,
        size=besthyper$size,
        decay=besthyper$decay,
        MaxNWts=Inf,entropy=T,maxit=1000,
        reltol = 1.0e-8,trace=F)
    
    predmat1   <- t(ktsppredmat(dat,lst.ktsp[[1]]$ktsp))
    predmat2   <- t(ktsppredmat(dat,lst.ktsp[[2]]$ktsp))
    predmat3   <- t(ktsppredmat(dat,lst.ktsp[[3]]$ktsp))
    predmat4   <- t(ktsppredmat(dat,lst.ktsp[[4]]$ktsp))
    predmat5   <- t(ktsppredmat(dat,lst.ktsp[[5]]$ktsp))
    
    pred.nn1 <- predict(nn.models[[1]],predmat1[, 1:5])
    pred.nn2 <- predict(nn.models[[2]],predmat2[, 1:5])
    pred.nn3 <- predict(nn.models[[3]],predmat3[, 1:5])
    pred.nn4 <- predict(nn.models[[4]],predmat4[, 1:5])
    pred.nn5 <- predict(nn.models[[5]],predmat5[, 1:5])
    
    preds <- (pred.nn1+pred.nn2+pred.nn3+pred.nn4+pred.nn5) >2.5
    auc.nn1 <- evalfun(preds>0.5,grp)
    auc.nn1 <- evalfun(pred.nn1>0.5,grp)
    auc.nn2 <- evalfun(pred.nn2>0.5,grp)
    auc.nn3 <- evalfun(pred.nn3>0.5,grp)
    auc.nn4 <- evalfun(pred.nn4>0.5,grp)
    auc.nn5 <- evalfun(pred.nn5>0.5,grp)
    
    pred.nn <- predict(nn.models[[1]],predmat)
    auc.nn <- evalfun(pred>0.5,grp)
    
    
    combineFunc <- function(x){predict(model,x[seq_len(besthyper$features)]) > 0.5}
    
    combineFunc <- function(x){predict(nn.models[[2]],x[seq_len(besthyper$features)]) > 0.5}
    pred.ktspnn <- predict(ktsp,dat,combineFunc)
    auc.ktspnn  <- evalfun(pred.ktspnn,grp)
    
    pred.ktspnn <- predict(ktsp,expr.test.t.knn,combineFunc)
    auc.ktspnn  <- evalfun(pred.ktspnn,clin.test$her2)
    
    
    return(list(hyperparams=pf,combineFunc=combineFunc,ktsp=ktsp,
                bestparams=besthyper))
}


trainbestnn <- function(besthyper,lst.ktsp,dat,grp){
    models = lapply(seq_along(lst.ktsp), function(i){
        x <- lst.ktsp[[i]]
        dat.train <- dat[ ,x$idx.train, drop=F]
        grp.train <- grp[x$idx.train]
        dat.test  <- dat[ ,x$idx.test, drop=F]
        grp.test  <- grp[x$idx.test]
        ktsp      <- x$ktsp 
        predmat.train   <- t(ktsppredmat(dat.train,ktsp))
        predmat.test   <- t(ktsppredmat(dat.test,ktsp))
        
        model <- nnet(
            predmat.train[ ,seq_len(besthyper$features),drop=F], 
            grp.train,
            size=besthyper$size,
            decay=besthyper$decay,
            MaxNWts=Inf,entropy=T,maxit=1000,
            reltol = 1.0e-8,trace=F)
    })
    return(models)
}

cvevalnn <- function(dat,grp,lst.ktsp,evalfun,size,features,decay){
    params = expand.grid(size,features,decay)
    colnames(params) <- c("size","features","decay")
    mat.out <- matrix(0,nrow=nrow(params)*length(lst.ktsp),ncol=6)
    #iterate over cv splits
    
    rowidx = 1
    
    #iter over cross validation splits
    for(i in seq_along(lst.ktsp)){
        x <- lst.ktsp[[i]]
        dat.train <- dat[ ,x$idx.train, drop=F]
        grp.train <- grp[x$idx.train]
        dat.test  <- dat[ ,x$idx.test, drop=F]
        grp.test  <- grp[x$idx.test]
        ktsp      <- x$ktsp 
        predmat.train   <- t(ktsppredmat(dat.train,ktsp))
        predmat.test   <- t(ktsppredmat(dat.test,ktsp))
        
        #iterate over hyperprams
        for(j in seq_len(nrow(params))){
            cat(".")
            model <- nnet(
                predmat.train[ ,seq_len(params$features[j]),drop=F], 
                grp.train,size=params$size[j],
                decay=params$decay[j],MaxNWts=Inf,entropy=T,maxit=1000,
                reltol = 1.0e-8,trace=F)
            pred.test <- predict(model,
                                 predmat.test[ ,seq_len(params$features[j]),drop=F])
            pred.train <- predict(model,
                                  predmat.train[ ,seq_len(params$features[j]),drop=F])
            
            pred.ptest <- predict(model,
                                  ptest[ ,seq_len(params$features[j]),drop=F])
            mat.out[rowidx,1]  <- evalfun(pred.train,grp.train)
            mat.out[rowidx,2]  <- evalfun(pred.test,grp.test)
            mat.out[rowidx,3]  <- i
            mat.out[rowidx,4]  <- params$size[j]
            mat.out[rowidx,5]  <- params$features[j]
            mat.out[rowidx,6]  <- params$decay[j]
            rowidx = rowidx + 1        
        }  
    }
    mat.out <- as.data.frame(mat.out)
    colnames(mat.out) <- c("Train","Test","cv","size","features","decay")
    #sum out cv splits
    pf <- mat.out %.% group_by(size,features,decay) %.%
        summarise(Train=mean(Train),Test=mean(Test))    
    return(pf)
}
