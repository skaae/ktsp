#' Shrink a KTSP to given size
#' @param ktsp a KTSP
#' @param dat training data
#' @param k either an integer specifying the new size of the KTSP or a vector 
#' specifying wich TSP's should be kept. 
#' @export
#' @examples
#' data(ktspdata)
#' ktsp <- ktspcalc(dat,grp,5,"My KTSP")
#' ktspplotpairs(ktsp,dat)
#' ktsp.small <- ktspshrink(ktsp,dat,c(1,3,5))
#' ktspplotpairs(ktsp.small,dat,file="example.pdf")
ktspshrink <- function(ktsp,dat,k){
    stopifnot(class(ktsp) =="ktsp")
        stopifnot(k <= ktsp$k)
        k.val = k
        k = seq_len(k)
    

    ktsp$TSPs <- ktsp$TSPs[k, , drop=F]
    ktsp$score <- ktsp$score[k,drop=F]
    ktsp$geneNames <- ktsp$geneNames[k, , drop=F]
    ktsp$k <- k.val
    ktsp$train.pred <- predict(ktsp,dat)
    return(ktsp)
}