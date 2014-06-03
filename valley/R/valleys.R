

#' Wrapper function arround valleysG to find valleys of an R object
#' @param x numeric vector of node intensities to consider.
#' @param edges integer matrix of 2 columns defining undirected edges of the graph. 
#'        If NULL, the edges net is generated from the structure of x depending if it is a vector or a matrix.
#' @return an integer vector of size length(x) pointing to the valley of each node
#' @export
#' @author Julien Prados
#' @examples
#'  set.seed(999)
#'  y <- rnorm(50)
#'  x <- seq_along(y)
#'  v <- valleys(y)
#'  v[is.na(v)] <- which.min(y)
#'  print(v)
#'  plot(x,y,type="o",pch="+")
#'  text(x,y,v,adj=c(0.5,1.5),cex=0.75)
#'  segments(x,y,x,y[v],col="blue",lwd=3)
#'  
#'  data(mir)
#'  y <- mir[seq(1,nrow(mir)-200,by=10),seq(1,ncol(mir),by=10),]
#'  y <- rowMeans(y,dims=2)
#'  plot(c(1,ncol(y)),c(1,nrow(y)),type="n")
#'  rasterImage(y,1,1,ncol(y),nrow(y),interpolate=FALSE)
#'  y <- y[nrow(y):1,]
#'  v <- valleys(-y)
#'  v[is.na(v)] <- which.min(-y)
#'  d <- (-y + y[v])
#'  msk <- d>0.1
#'  points(col(y)[msk],row(y)[msk],pch=2,col="blue")
#'  text(col(y)[msk],row(y)[msk],format(d[msk],digits=2),adj=c(0.5,-0.5),cex=0.75)
valleys <- function(x,edges=NULL) {
  if(is.null(edges)) {
    if (is.matrix(x)) {
      i <- matrix(seq_along(x),nrow(x),ncol(x))
      edges <- cbind(c(i[-nrow(i),],i[,-ncol(i)]),c(i[-1L,],i[,-1L]))
    } else {
      to <- seq_along(x)[-1]
      edges <- cbind(to-1,to)
    }
  }
  valleysG(as.numeric(x),edges)
}


