#'  tiedRank(x) applies rank to the columns of the matrix x. 
#' 
#'  tiedRank(x) applies rank to the columns of the matrix x. It is used to calculate the secondary score.
#'
#' @param x a data matrix
#'
#' @return x a rank matrix
#'
#' @keywords rank
#'
tiedRank<-function(x)
{
    return((apply(x,2,rank)));

}