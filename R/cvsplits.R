#' Create Crossvalidation split indices
#' 
#' @param Nobs Number of observations in dataset
#' @param K fold cross validation
#' @param seed seed for cross validation splits
#' @return a list of length ??code(k) with train and test indices
#' @export
cvsplits <- function(Nobs,K=5,seed=1234){
    set.seed(seed)
    rs <- runif(Nobs)
    id <- seq(Nobs)[order(rs)]
    k <- as.integer(Nobs*seq(1,K-1)/K)
    k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
    k[,1] <- k[,1]+1
    l <- lapply(seq.int(K),function(x,k,d) 
        list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
             test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
    return(l)   
}