#' Calculate prediction from KTSP object
#' 
#' Takes a KTSP and data as input and returns the predictions as a vector.
#' Optionally takes a function for combining output from the k-tsp's into a 
#' single prediction for each sample. If no combine function is given 
#' majority prediction is used.
#' 
#' @param object A KTSP object
#' @param data a matrix with rows as genes and samples as columns
#' @param combineFunc A function for combining the k-tsp's into a single 
#' prediction. Default is majority predict. See \code{\link{majorityPredict}} 
#' for example
#' @param ... not used
#' @rdname predict
#' @return vector of predictions
#' @S3method predict ktsp
#' @examples
#' data(ktspdata)
#' ktsp <- ktspcalc(dat,grp,5)
#' pred <- predict(ktsp,dat)
#' @export
predict.ktsp<-function(object, data, combineFunc=majorityPredict,...) {
    stopifnot(is.function(combineFunc))
    ktsp <- object
    tspIndexes <- t(apply(ktsp$geneNames, 1, FUN <- function(x,y = data) {
        a <- which(rownames(y) %in% x[1])
        b <- which(rownames(y) %in% x[2])
        c(a, b)
    }))
    if (is.vector(data)) {
        x <- data[tspIndexes[, 1]] < data[tspIndexes[, 2]]
        s <- combineFunc(x)
    }
    else {
        x <- data[tspIndexes[, 1], ] < data[tspIndexes[, 2], 
                                            ]
        if(is.null(dim(x))){
            s <- sapply(x,combineFunc)
        }else{
            s <- apply(as.matrix(x), 2, combineFunc)
        }
    }
    return(as.integer(s))
}