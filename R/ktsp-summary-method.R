#' Print summary statistics
#' @method summary ktsp
#' @import caret 
#' @param object ktsp object
#' @param ... not used
#' @S3method summary ktsp
summary.ktsp <- function(object,...){
    ktsp <- object
    pred <- ktsp$train.pred
    grp <- ktsp$train.grp
    print(ktsp)
    cat("---------------------------------\n")
    cat(paste0("AUC: ",CalcAUC(pred,grp),"\n"))
    print(confusionMatrix(as.factor(pred),reference=as.factor(grp)))
    cat("---------------------------------\n")
}