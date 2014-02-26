require(e1071)
ktsptrainsvm <- function(dat, grp, k.max, cross, evalfun = CalcAUC,
                        splits = NULL, seed = 1234,c,gamma,features,precalc=NULL){
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
    
    pf       <- cvevalsvm(dat,grp,lst.ktsp,evalfun,size=size,features=features,
                         decay=decay)
    besthyper <- pf[which.max(pf$Test), ,drop=F]
    
    ktsp      <- ktspcalc(dat,grp,besthyper$features,"test")
    predmat   <- t(ktsppredmat(dat,ktsp))
    svm    <- svm(x=predmat,
                  y=grp,
                  gamma=besthyper$gamma,
                  cost=besthyper$c,
                  type="C-classification")
    combineFunc <- function(x){predict(svm,x[seq_len(besthyper$features)])}
    
    return(list(hyperparams=pf,combineFunc=combineFunc,ktsp=ktsp,
                bestparams=besthyper))
}


cvevalsvm <- function(dat,grp,lst.ktsp,evalfun,c,gamma,features){
    params = expand.grid(c,gamma,features)
    colnames(params) <- c("c","gamma","features")
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
            model <- svm(
                predmat.train[ ,seq_len(params$features[j]),drop=F], 
                grp.train,type="C-classification",
                gamma=params$gamma[j],
                cost=params$c[j])
            pred.test <- predict(model,
                                 predmat.test[ ,seq_len(params$features[j]),drop=F])
            pred.train <- predict(model,
                                  predmat.train[ ,seq_len(params$features[j]),drop=F])
            mat.out[rowidx,1]  <- evalfun(pred.train,grp.train)
            mat.out[rowidx,2]  <- evalfun(pred.test,grp.test)
            mat.out[rowidx,3]  <- i
            mat.out[rowidx,4]  <- params$c[j]
            mat.out[rowidx,5]  <- params$gamma[j]
            mat.out[rowidx,6]  <- params$features[j]
            rowidx = rowidx + 1        
        }  
    }
    mat.out <- as.data.frame(mat.out)
    colnames(mat.out) <- c("Train","Test","cv","c","gamma","features")
    #sum out cv splits
    pf <- mat.out %.% group_by(c,gamma,features) %.%
        summarise(Train=mean(Train),Test=mean(Test))
    
    return(pf)
}
