

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
#'  plot(c(1,ncol(mir)),c(1,nrow(mir)),type="n")
#'  rasterImage(mir,1,1,ncol(mir),nrow(mir),interpolate=FALSE)
#'  mir <- mir[nrow(mir):1,]
#'  v <- valleys(-mir)
#'  v[is.na(v)] <- which.min(-mir)
#'  d <- (-mir + mir[v])
#'  msk <- d>0.1
#'  points(col(mir)[msk],row(mir)[msk],pch=2,col="blue")
#'  text(col(mir)[msk],row(mir)[msk],format(d[msk],digits=2),adj=c(0.5,-0.5),cex=0.75)
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


