#' Calculate prediction from ktsp
#' 
#' Returns a matrix of size samples by number_of_tsp
#' Each row is the TSP prediction of for a sample
#' 
#' @param data expresssion matrix. Genes are rows and samples are columns
#' @param ktsp ktsp
#' @export
#' @return matrix of predictions. size samples by number_of_tsp
#' @examples
#' #Example 1
#' data(ktspdata)
#' ktsp <- ktspcalc(dat,grp,5,name="Test TSP")
#' predmat <- ktsppredmat(dat,ktsp)
#' pred <- apply(predmat,2,majorityPredict)
#' AUC <- CalcAUC(pred,ktsp$train.grp)
#' stopifnot(AUC==CalcAUC(ktsp$train.pred,ktsp$train.grp))
#' 
#' #Example 2 - load data and precalculated tsp
#' data(her2) 
#' data(ktsp)
#' predmat <- ktsppredmat(expr.train.t.knn,ktsp)
ktsppredmat <- function(data, ktsp){
    tspIndexes <- t(apply(ktsp$geneNames, 1, FUN <-function(x,y = data) {
        a <- which(rownames(y) %in% x[1])
        b <- which(rownames(y) %in% x[2])
        c(a, b)
    }))
    x <- data[tspIndexes[, 1], ] < data[tspIndexes[, 2], ]
    rownames(x) <- paste0("TSP",seq(1,ktsp$k))
    x <- x*1
    return(x)
}