#' Print KTSP
#' 
#' Print basic decription of KTSP class
#' @S3method print ktsp
#' @param x a ktsp
#' @param y unused
#' @param ... unused
#' @examples
#' data(ktspdata)
#' ktsp <- ktspcalc(dat,grp,5,"Estrogen Receptor")
#' print(ktsp)
print.ktsp <- function(x,y,...){
    cat(paste0("A KTSP object with ",x$k,"-TSP's. Name of TSP: ",x$name))
}