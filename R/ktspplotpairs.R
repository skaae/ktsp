#' Plot KTSP pairs
#' 
#' Plots all pairs in a TSP. Use \code{\link{ktspshrink}} if you do only 
#' want to plot some of the pairs 
#' @param ktsp a KTSP
#' @param dat training data
#' @param file optionally set save path. 
#' @export
#' @import ggplot2
#' @examples
#' data(ktspdata)
#' ktsp <- ktspcalc(dat,grp,5,"My KTSP")
#' ktspplotpairs(ktsp,dat)
#' ktspplotpairs(ktsp,dat,file="example.pdf")
ktspplotpairs <- function(ktsp,dat,file=NULL){
    stopifnot(class(ktsp) == "ktsp")
    stopifnot(ncol(dat) == length(ktsp$train.grp))
    name <- ktsp$name
    if(name==""){
        name=paste0("Expression of genes in ",ktsp$k,"-TSP")
    }
    
    l = nrow=ncol(dat)
    pdat = as.data.frame(matrix(0,l*ktsp$k,ncol=6))
    colnames(pdat) <- c("x","y","Label","k","Gene1","Gene2")
    for(i in seq(1,ktsp$k)){
        r = seq(((i-1)*l+1),i*l)
        pdat[r,1] = dat[ktsp$TSPs[i,1], ]
        pdat[r,2] = dat[ktsp$TSPs[i,2], ]
        pdat[r,3] = ktsp$train.grp
        pdat[r,4] = paste("TSP",i)
        pdat[r,5] = paste(paste("TSP",i,"-"),"X:",ktsp$geneNames[i,1])
        pdat[r,6] = paste("Y:",ktsp$geneNames[i,2])
    }
    pdat$Label <- as.factor(pdat$Label)
    
    
    p <- ggplot(pdat,aes(x=x,y=y,color=Label)) +
        geom_point() +
        facet_wrap(Gene1 ~ Gene2,ncol=3) +
        scale_color_brewer(palette="Set1",type="qual") +
        theme(strip.text.x = element_text(size = rel(0.7))) +
        ggtitle(name) +
        ylab("Expression second gene") +
        xlab("Expression first gene")
    
    if(is.null(file)){
        print(p)
    }else{
        print(p)
        pdf(file)
        print(p)
        dev.off()
    }
    
    
}
