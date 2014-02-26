#' Create RandomForest combine function
#' 
#' Using CV this function trains a RandomForest ontop of a KTSP. User specify 
#' hyper parameter space to be searched.
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param k.max largest TSP to evaluate
#' @param cross number of crosss validation folds 
#' @param evalfun function for eval performance on CV splits. Defaults to AUC 
#' @param splits optionally supply the cv split indices
#' @param seed seed for cross validation splits
#' @param mtry randomforest hyperparam \code{link{randomForest}}. Specify values
#' to be tried as vector
#' @param nodesize randomfforest hyperparam \code{link{randomForest}}. Specify 
#' values to be tried as vector
#' @param features vector of number of features to try
#' see \code{\link{CalcAUC}}
#' @export
#' @import dplyr
#' @import randomForest
#' @return A list with the following fields
#'      \item{hyperparams}{A matrix with average test/trainset performance using
#'      different hyperparameter values}
#'      \item{bestparams}{best hyper parameters }
#'      \item{ktsp}{A ktsp trained on the training data using best k tsp's}
#'      \item{combineFunc}{The combine function used/created}
#' @family hyperparams
#' @examples
#' data(ktspdata)
#' cross = 5
#' splits <- cvsplits(ncol(dat), cross)
#' res <- ktsptrainrf(dat,grp,100,cross,splits=splits, mtry=c(1,2,3),
#'                    nodesize=c(1,2,3,4,5),features=c(10,50,100))
#' pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
#' summary(besttsp)
#' 
#' data(her2)
#' cross = 5
#' k.max = 200
#' res <- ktsptrainrf(expr.train.t.knn,clin$her2,k.max,cross,splits=splits, 
#'                    mtry=c(1,2,3),nodesize=c(1,2,3,4,5),
#'                    features=c(5,10,30,50,100,200))
#' pred <- predict(res$ktsp,expr.train.t.knn,combineFunc=res$combineFunc)
#' auc <- CalcAUC(res$ktsp,clin.test$her2)
ktsptrainrf <- function(dat, grp, k.max, cross, evalfun = CalcAUC,
                        splits = NULL, seed = 1234,
                        mtry=NULL,nodesize=NULL,features=NULL,precalc=NULL){
    stopifnot(all(!c(is.null(mtry),is.null(nodesize),is.null(features))==T))
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
    pf       <- cvevalrf(dat,grp,lst.ktsp,mtry,nodesize,features,evalfun)
    besthyper <- pf[which.max(pf$Test), ,drop=F]
    
    ktsp      <- ktspcalc(dat,grp,besthyper$features,"test")
    predmat   <- t(ktsppredmat(dat,ktsp))
    rf    <- randomForest(x=predmat,
                          y=factor(grp),ntree=500,mtry=besthyper$mtry,
                          nodesize=besthyper$nodesize)
    combineFunc <- function(x){predict(rf,x[seq_len(besthyper$features)])}
    
    return(list(hyperparams=pf,combineFunc=combineFunc,ktsp=ktsp,
                bestparams=besthyper))
}


#' Evaluate random forest performance on cv splits
#' 
#' @param lst.ktsp list of ktsp, format as output of \code{\link{cvtrain}}
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param evalfun evaluation function see \code{\link{CalcAUC}}
#' @param mtry mtry hyperparam vals
#' @param nodesize nodesize hyper param vals
#' @param features hyperparam vals
#' @return df with the performance for TSP's up to k. Columns:
#' #'   \item{Train}{Training performance}
#'      \item{Test}{Test performance}
#'      \item{k}{Number of included tsp's}
cvevalrf <- function(dat,grp,lst.ktsp,mtry,nodesize,features,evalfun){
    params = expand.grid(mtry,nodesize,features)
    colnames(params) <- c("mtry","nodesize","features")
    
    #remove hyperparams where features > mtry
    params[params$mtry > params$features,"mtry"] <- 
                               params[params$mtry > params$features,"features"]
    params = as.data.frame(params[!duplicated(params), , drop=F])
    
    mat.out <- matrix(0,nrow=nrow(params)*length(lst.ktsp),ncol=6)
    #iterate over cv splits
    
    rowidx = 1
    
    #iter over cross validation splits
    for(i in seq_along(lst.ktsp)){
        x <- lst.ktsp[[i]]
        dat.train <- dat[ ,x$idx.train, drop=F]
        grp.train <- factor(grp[x$idx.train])
        dat.test  <- dat[ ,x$idx.test, drop=F]
        grp.test  <- factor(grp[x$idx.test]) 
        ktsp      <- x$ktsp 
        predmat.train   <- t(ktsppredmat(dat.train,ktsp))
        predmat.test   <- t(ktsppredmat(dat.test,ktsp))
        #iterate over hyperprams
        for(j in seq_len(nrow(params))){
            cat(".")
            forest <- randomForest(
                predmat.train[ ,seq_len(params$features[j]),drop=F], 
                grp.train,ntree=500,mtry=params$mtry[j],
                nodesize=params$nodesize[j])
            pred.test <- predict(forest,
                                 predmat.test[ ,seq_len(params$features[j]),drop=F])
            pred.train <- predict(forest,
                                  predmat.train[ ,seq_len(params$features[j]),drop=F])
            mat.out[rowidx,1]  <- evalfun(pred.train,grp.train)
            mat.out[rowidx,2]  <- evalfun(pred.test,grp.test)
            mat.out[rowidx,3]  <- i
            mat.out[rowidx,4]  <- params$mtry[j]
            mat.out[rowidx,5]  <- params$nodesize[j]
            mat.out[rowidx,6]  <- params$features[j]
            rowidx = rowidx + 1        
        }  
    }
    mat.out <- as.data.frame(mat.out)
    colnames(mat.out) <- c("Train","Test","cv","mtry","nodesize","features")
    #sum out cv splits
    pf <- mat.out %.% group_by(mtry,nodesize,features) %.%
        summarise(Train=mean(Train),Test=mean(Test))
    
    return(pf)
}







