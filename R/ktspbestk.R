#' Determine best K using cross validation
#' 
#' Calculates the optimal number of TSP's using cross validation
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels, vector
#' @param k.max largest KTSP to evaluate, may be omitted if precalc is supplied
#' @param cross number of crosss validation folds,may be omitted if precalc 
#' is supplied
#' @param file optional file name for plot. Ignored if fs =TRUE
#' @param seed seed for cross validation splits
#' @param evalfun function for evaluating performance on CV splits. Defaults to AUC 
#' @param name Optional name for plot 
#' @param splits optionally supply the cv split indices, if supplied cross and
#' seed is ignored
#' @param type Feature selection type. If "seq" the function
#' tries TSP's sequentially, i.e starting with a ktsp with TSP-1, then a KTSP
#' with TSP-1+2, 1+2+3 etc. If "forward" greedy forward search, if "backward", 
#' greedy backward search
#' @param precalc optionally supply a precalculated list of CV tsps, output from
#' \code{\link{cvtrain}}, to reduce computational time. If supplied k.max,
#' cross and seed are ignored.
#' @param reg regularization function. see \code{\link{Reg}}
#' @export
#' @import dplyr
#' @import ggplot2
#' @return A list with the following fields
#'      \item{features}{either bestk value (fs = "seq") or the selcted features
#'                     (fs = "backward/forward")}
#'      \item{ktsp}{A ktsp trained on the training data using best k tsp's}
#'      \item{combineFunc}{The combine function used/created}
#' @family hyperparams
#' @examples
#' \dontrun{
#' # Example 1 - selecting top TSP's (no FS)
#' data(ktspdata)
#' res <- ktspbestk(dat,grp,20,3,file="example.pdf",name="Test TSP",
#'                 type="seq")
#' pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
#' summary(res$ktsp)
#' 
#' # forward search
#' res <- ktspbestk(dat,grp,20,3,file="example.pdf",name="Test TSP",
#'                 type="forward")
#' pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
#' 
#' # backward search
#' res <- ktspbestk(dat,grp,20,3,file="example.pdf",name="Test TSP",
#'                 type="backward")
#' pred <- predict(res$ktsp,dat,combineFunc=res$combineFunc)
#' 
#' # Example - her2
#' require(datathesis)
#' data(her2) # load expression/clinical data and precalculated ktsp list
#' resFSseq<- ktspbestk(expr.train.t.knn,clin.train$her2,type="seq",
#'                        file="example.pdf",name="HER2 FS:sequential",
#'                        precalc=lst.ktsp)             
#' resFSf <- ktspbestk(expr.train.t.knn,clin.train$her2,type="forward",
#'                        file="example.pdf",name="HER2 FS:Forward",
#'                        precalc=lst.ktsp)  
#' resFSb <- ktspbestk(expr.train.t.knn,clin.train$her2,type="backward",
#'                        file="example.pdf",name="HER2 FS:Backward",
#'                        precalc=lst.ktsp)  
#'                                     
#' predFSseq  <- predict(resFSseq$ktsp,expr.test.t.knn,
#'                         combineFunc=resFSseq$combineFunc)
#' predFSf  <- predict(resFSf$ktsp,expr.test.t.knn,
#'                          combineFunc=resFSf$combineFunc)
#' predFSf  <- predict(resFSb$ktsp,expr.test.t.knn,
#'                          combineFunc=resFSb$combineFunc)
#' aucFSseq <- CalcAUC(predFSseq,clin.test$her2)
#' aucFSf  <- CalcAUC(predFSf,clin.test$her2)
#' aucFSb  <- CalcAUC(predFSb,clin.test$her2)
#' }
ktspbestk <-function(dat,grp,k.max,cross,type,evalfun=CalcAUC,name=NULL,splits=NULL,
                     file=NULL,seed=1234,precalc=NULL,reg=noReg){    
    stopifnot(length(grp) == ncol(dat))
    # train a ktsp on each of the cross validation splits
    if(is.null(k.max)){
        k.max <- lst.ktsp[[1]]$ktsp$k
    }
    
    
    if(is.null(precalc)){
        if(is.null(splits)){
            splits   <- cvsplits(ncol(dat),cross,seed)
        }
        lst.ktsp <- cvtrain(splits,dat,grp,k.max,name)
    }else{
        lst.ktsp <- precalc
        cross = length(lst.ktsp)
    }  
    if(type=="seq"){
        pf       <- cvevalbestk(dat,grp,lst.ktsp,evalfun=evalfun,reg=reg,k.max=k.max)
        pf <- pf[1:k.max, ,drop=F]
        pdat     <- cvcreatepdatbestk(pf)
        cvplotbestk(pdat,name,file)
        bestk <- pf$k[which.max(pf$Test)]
        cat("train best KTSP\n")
        ktsp <- ktspcalc(dat,grp,bestk)
        cat("-----------------Selected Genes------------------\n")
        print(ktsp$geneNames)
        cat("-------------------------------------------------\n")
        return(list(features=bestk,combineFunc=majorityPredict,ktsp=ktsp))
    }else{
        stop()
        stopifnot(type %in% c("forward","backward"))
        temp <- cvevalbestkfs(dat,grp,lst.ktsp,evalfun=evalfun,type=type)
        cat("-----------------Selected TSP's------------------\n")
        print(temp$features)
        cat("-------------------------------------------------\n")
        cat("train best KTSP\n")
        ktsp <- ktspcalc(dat,grp,temp$k)
        idx <- sapply(temp$features,function(x) as.numeric(gsub("TSP","",x)))
        feat <- ktsp$geneNames[idx, , drop=F]
        return(list(features=feat,combineFunc=temp$combineFunc,ktsp=ktsp))
        
    }
    
    
}

#' Evaluate performance seq search CV
#' 
#' @param lst.ktsp list of ktsp, format as output of \code{\link{cvtrain}}
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param evalfun evaluation function see \code{\link{CalcAUC}} for example
#' @return df with the performance for TSP's up to k. Columns:
#'  \item{Train}{Training performance}
#'  \item{Test}{Test performance}
#'  \item{k}{Number of included tsp's}
cvevalbestk <- function(dat,grp,lst.ktsp,evalfun,reg,k.max){
    cat("Feature selection: Try TSP's sequentially \n")
    k.max <- lst.ktsp[[1]]$ktsp$k
    k.seq <- seq(1,k.max)
    stopifnot(sort(unique(grp)) == c(0,1))
    mat.out <- matrix(0,nrow=length(k.seq)*length(lst.ktsp),ncol=3)
    rowidx = 1
    for(i in seq_along(lst.ktsp)){  # cross validation splts
        x <- lst.ktsp[[i]]
        dat.train <- dat[ ,x$idx.train, drop=F]
        grp.train <- grp[x$idx.train] 
        dat.test  <- dat[ ,x$idx.test, drop=F]
        grp.test  <- grp[x$idx.test] 
        ktsp      <- x$ktsp  
        
        for(k.new in k.seq){
            ktsp.new           <- ktspshrink(ktsp,dat.train,k.new)
            pred.train         <- predict(ktsp.new,dat.train) 
            pred.test          <- predict(ktsp.new,dat.test)
            regCost            <- reg(k.new,k.max)
            mat.out[rowidx,1]  <- evalfun(pred.train,grp.train) -regCost
            mat.out[rowidx,2]  <- evalfun(pred.test,grp.test)   -regCost
            mat.out[rowidx,3]  <- k.new
            rowidx = rowidx + 1        
        }  
    }
    mat.out <- as.data.frame(mat.out)
    colnames(mat.out) <- c("Train","Test","k")
    pf <- mat.out %.% group_by(k) %.%
        summarise(Train=mean(Train),Test=mean(Test)) %.%
        arrange(k)
    
    return(pf)
}

#' Select bestk using forward/backward feature selection
#' 
#' @param lst.ktsp list of ktsp, format as output of \code{\link{cvtrain}}
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param evalfun evaluation function see \code{\link{CalcAUC}} for example
#' @param type "forward" or "backward" search
#' @return list with the following fields
#'  \item{combineFunc}{CombineFunc to be used with a ktsp trained with k tsp's}
#'  \item{k}{train ktsp with the many features}
#'  \item{features}{selected features}
# cvevalbestkfs <- function(dat,grp,lst.ktsp,evalfun,type){
#     if(type=="forward"){
#         cat("Feature selection: forward search \n")
#         featureselect <- forward.search
#     }else{
#         cat("Feature selection: backward search \n")
#         featureselect <- backward.search
#     }
#     cross = length(lst.ktsp)
#     
#     # precalc prediction matrices for test data
#     predmats <- lapply(lst.ktsp,function(x){
#         ktsp <- x$ktsp
#         dat.test  <- dat[ ,x$idx.test, drop=F]
#         predmat <- ktsppredmat(dat.test,ktsp)
#         grp     <- grp[x$idx.test] 
#         return(list(predmat=predmat,grp=grp))
#     })
#     
#     #create evaluation function
#     this.eval <- function(selfeat) evalfs(selfeat,cross,dat,grp,
#                                           predmats,evalfun)
#     ntsp <- nrow(lst.ktsp[[1]]$ktsp$TSPs)
#     attrib <- paste0("TSP",seq_len(ntsp))
#     
#     temp <- featureselect(attrib,this.eval)
#     idx <- sapply(temp$features,function(x) as.numeric(gsub("TSP","",x)))
#     combineFunc <- function(x) {
#        
#         x <- x[idx]
#         pred <- mean(x) > 0.5
#         return(pred)
#     }
#     return(list(combineFunc=combineFunc, k=ntsp,features=features))
# }




# cvevalbestkfs <- function(dat,grp,lst.ktsp,evalfun,type){
#     if(type=="forward"){
#         cat("Feature selection: forward search \n")
#         featureselect <- forward.search
#     }else{
#         cat("Feature selection: backward search \n")
#         featureselect <- backward.search
#     }
#     cross = length(lst.ktsp)
#     ntsp <- nrow(lst.ktsp[[1]]$ktsp$TSPs)
#     attrib <- paste0("TSP",seq_len(ntsp))
#     
#     mat.out <- matrix(0,nrow=length(k.seq)*length(lst.ktsp),ncol=3)
#     rowidx = 1
#     for(i in seq_along(lst.ktsp)){  # cross validation splts
#         x <- lst.ktsp[[i]]
#         dat.train <- dat[ ,x$idx.train, drop=F]
#         grp.train <- grp[x$idx.train] 
#         dat.test  <- dat[ ,x$idx.test, drop=F]
#         grp.test  <- grp[x$idx.test] 
#         ktsp      <- x$ktsp 
#         
#         predmat.test <- ktsppredmat(dat.test,ktsp)
#         predmat.train <- ktsppredmat(dat.train,ktsp)
# 
#         this.eval <- function(selfeat) evalfs(selfeat,predmat.test,
#                                               grp.test,evalfun)
#         temp <- featureselect(attrib,this.eval)
#         print(temp$features)
#     }
# #         for(k.new in k.seq){
# #             ktsp.new   <- ktspshrink(ktsp,dat.train,k.new)
# #             pred.train <- predict(ktsp.new,dat.train) 
# #             pred.test          <- predict(ktsp.new,dat.test)
# #             mat.out[rowidx,1]  <- evalfun(pred.train,grp.train)
# #             mat.out[rowidx,2]  <- evalfun(pred.test,grp.test)
# #             mat.out[rowidx,3]  <- k.new
# #             rowidx = rowidx + 1        
# #         }  
#     }
#     mat.out <- as.data.frame(mat.out)
#     colnames(mat.out) <- c("Train","Test","k")
#     pf <- mat.out %.% group_by(k) %.%
#         summarise(Train=mean(Train),Test=mean(Test)) %.%
#         arrange(k)
#     
#     return(pf)
# }

#' evalute performance, used in forward/backward search
#' 
#' @param predmats list of precalculated prediction matrices
#' @param dat gene expression df Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param cross number of cross validation splits
#' @param attrib features to try
#' @param evalfun evaluation function see \code{\link{CalcAUC}} for example
#' @return list with the following fields
#'  \item{combineFunc}{CombineFunc to be used with a ktsp trained with k tsp's}
#'  \item{k}{train ktsp with the many features}
#'  \item{features}{selected features}
evalfs <- function(attrib,predmat,grp,evalfun){
    predmat <- predmat[attrib, , drop=F]
    pred <- 1*apply(predmat,2,majorityPredict) # get majority predict
    perf <- evalfun(pred,grp)
    return(perf)             # return performance
}

#' Create plotting data
#' 
#' @param pf output from  \code{\link{cvevalbestk}}
#' @return a df in format suitable for ggplot2 plotting
cvcreatepdatbestk <- function(pf){
    
    bestk.idx <- which.max(pf$Test)
    pdat.max=data.frame(
        y=c(pf$Test[bestk.idx]),
        x=c(pf$k[bestk.idx]),
        value=c(paste0(round(pf$Test[bestk.idx],2),
                       ", k=",pf$k[bestk.idx])))
    pdat <- melt(pf,id.vars=c("k"))
    colnames(pdat) <- c("k","Dataset","value")
    return(list(pdat=pdat,pdat.max=pdat.max))
}

#' Plot cv performance on training and test set
#' 
#' @param pdat output from  \code{\link{cvcreatepdatbestk}}
#' @param title plot title
#' @param file output file. If NULL no save
cvplotbestk <- function(pdat,title,file){
    pdat.max = pdat$pdat.max
    pdat = pdat$pdat
    p <- ggplot() +
        geom_line(data=pdat,aes(x=k,y=value,color=Dataset)) +
        geom_point(data=pdat,aes(x=k,y=value,color=Dataset)) +
        geom_point(data=pdat.max,aes(x=x,y=y)) +
        geom_text(data=pdat.max,aes(label=value,x=x,y=y),hjust=-0.1, vjust=0) +
        ggtitle(paste0(title,": TSP CV performance")) +
        xlab("Number of TSP's") +
        ylab("Performance") +
        scale_color_brewer(type="qual",palette="Set2",labels=c("Train","Test")) +
        theme(legend.position="bottom")
    
    if(is.null(file)){
        print(p)
    }else{
        print(p)
        pdf(file)
        print(p)
        dev.off()     
    }  
}




