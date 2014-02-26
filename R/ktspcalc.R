#' train KTSP classifier. 
#'
#' Train a KTSP classifier with n TSP's. 
#'
#' @param data gene matrix. Rows are genes, columns are samples
#' @param grp a vector of training labels encoded as 1 and 0
#' @param n The number of TSP's to include.
#' @param name optional name of KTSP. For plots etc
#' @return Retuns a KTSP object with the following fields
#'      \item{TSPs}{Matrix of indices of selected pairs. k by 2}
#'      \item{k}{Number of TSP's}
#'      \item{geneNames}{Matrix of names of selected pairs. k by 2}
#'      \item{name}{Optional name of ktsp used for plotting etc}
#'      \item{train.grp}{labels of training data}
#'      \item{train.pred}{prediction on training data}
#' @keywords train ktsp
#' @examples
#' data(ktspdata)
#' ktsp <- ktspcalc(dat,grp,5,"My KTSP")
#' summary(ktsp)
#' @export
ktspcalc<-function(data, grp, n,name=NULL){
    if(is.null(name)){
        name = paste0(n,"-TSP")
    }
        
    
    data = tiedRank( data );
    KTSPout =calculateSignedScore(  grp , data, data );
    TSPs = matrix(0,n,2);
    TSPscore = vector(length=n);
    TSPGenes = matrix("",n,2);
    geneNames = rownames(data);
    score = KTSPout$score;
    
    for(i in 1:n)
    {
        #find max of score matrix
        TSPs[i,] = arrayInd(which.max(score),.dim=dim(score));
        
        
        TSPscore[i] = score[TSPs[i,1],TSPs[i,2]]/2;
        TSPGenes[i,] = geneNames[TSPs[i,]];
        score[TSPs[i,1],]=0;
        score[TSPs[i,2],]=0;
        score[,TSPs[i,1]]=0;
        score[,TSPs[i,2]]=0;
    }
    ktsp <-list(TSPs=TSPs,
                 geneNames=TSPGenes,
                 k=n,
                 name=name,
                 train.grp = grp);
    class(ktsp) <- "ktsp"
    pred <- predict(ktsp,data)
    ktsp$train.pred <- pred
    
    return(ktsp)
}