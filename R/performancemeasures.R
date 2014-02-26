#' Calculate AUC score.
#' @param pred prediction vector
#' @param grp label vector
#' @return AUC
#' @import ROCR
#' @export
#' @family performance
CalcAUC <- function(pred,grp){
    require("ROCR")
    p   <- prediction(as.numeric(pred),grp)
    auc <- performance(p,"auc")
    auc <- unlist(slot(auc, "y.values"))
    return(auc)
}
#' Calculate Accuracy.
#' @param pred prediction vector
#' @param grp label vector
#' @return AUC
#' @import ROCR
#' @export
#' @family performance
CalcAcc <- function(pred,grp){
    acc = sum(pred==grp) / length(grp)
}
