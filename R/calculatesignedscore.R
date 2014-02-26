#' Calculates pairwise score between every gene in data1 and every gene in data2
#' 
#' Calculates pairwise score between every gene in data1 and every gene in data2. It is signed because
#' score(i,j)=P(X_i<X_j|1)-P(X_i<X_j|0) - P(X_i>X_j|1) +P(X_i>X_j|0) + C
#' The third and the forth terms are for avoiding the equlities. The C term is proportion to the secondary score
#' to avoid the ties.  C = (E(X_j-X_i|1)-E(X_j-X_i|0))/K where K is big enough to make sure that the secondary score
#' does not intervene with the primary score.
#' a,M,k,N are usable for fisher exact test. M is the size of the class situation == 0 and N is the total number
#' of samples. a_ij = #(X_i<X_j|0) k_ij = #(X_i<X_j)
#' P_ij = P(X_i<X_j|0) and Q(X_i<X_j|1)
#'
#' @param label group labels
#' @param data1 first gene matrix
#' @param data2 second gene matrix
#'
#' @return signedscore Signed socres
#'
#' @keywords Calculate signed score
#'
calculateSignedScore<-function(  label , data1, data2 )
{
    
    n = length(label);
    m1 = nrow(data1);
    m2 = nrow(data2);
    
    d<-
    .C(
    "CalculateSignedScoreCore",
    as.integer(label), as.integer(n),
    as.double(data1),as.integer(m1),as.double(data2),as.integer(m2),
    as.double(matrix(0,m1,m2)),as.double(matrix(0,m1,m2)),
    as.double(1),as.double(matrix(0,m1,m2)),
    as.double(1),as.double(matrix(0,m1,m2)),
    as.double(matrix(0,m1,m2))
    );
    
    retVal<-list(score=(matrix(d[[7]],nrow=m1)),a=(matrix(d[[7]],nrow=m1)),
                 M=(matrix(d[[8]],nrow=1)),k=(matrix(d[[9]],nrow=m1)),N=(matrix(d[[10]],nrow=1)),
                 P=(matrix(d[[11]],nrow=m1)),Q=(matrix(d[[12]],nrow=m1)));  

}