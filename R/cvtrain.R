#' Train a ktsp on each cv split
#' 
#' This function is used to evaluate CV performance e.g. in 
#' \code{\link{ktspbestk}} and \code{\link{ktsptrainrf}}
#' @param splits cv splits see \code{\link{cvsplits}}
#' @param dat gene expression matrix. Genes as rows and samples as columns
#' @param grp binary encoded labels
#' @param k number of tsp's
#' @param name TSP name
#' @export
#' @return list of length "number of corssval folds" with the following fields
#'      \item{ktsp}{a ktsp object}
#'      \item{idx.train}{indices of cv train data}
#'      \item{idx.test}{indices of cv test data}
cvtrain <- function(splits,dat,grp,k,name){
    lst.ktsp <- lapply(splits,function(spl){
        cat(sprintf("."))
        dat.train <- dat[ , spl$train, drop=F]
        grp.train <- grp[spl$train] 
        ktsp <- ktspcalc(dat.train,grp.train,k,name)
        return(list(ktsp=ktsp, idx.train=spl$train, idx.test=spl$test))
    })
    return(lst.ktsp)
}