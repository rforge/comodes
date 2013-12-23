setMethod(
  f="plot",
  signature = c("CoModes_res","numeric"),
  definition = function(x,y,...) {
    if (is.numeric(y)==FALSE)
      y <- 1:2
    
    if (length(y)!=2)
      y=1:2
    
    res.mca <- MCA(x@data@data,graph=FALSE)
    xlab <- paste("Axe",y[1],"of the multiple correspondance analysis")
    ylab <- paste("Axe",y[2],"of the multiple correspondance analysis")
    plot(res.mca$ind$coord[,y]+matrix(runif(2*nrow(res.mca$ind$coord[,y]),0,0.1),ncol=2),xlab=xlab,ylab=ylab,col=as.numeric(x@indiv@partition),pch=as.numeric(x@indiv@partition))
  }
)
# 
# setMethod(
#   f="plot",
#   signature = c("CoModes_res"),
#   definition = function(x...) {
#     y <- 1:2
#     res.mca <- MCA(x@data@data,graph=FALSE)
#     xlab <- paste("Axe",y[1],"of the multiple correspondance analysis")
#     ylab <- paste("Axe",y[2],"of the multiple correspondance analysis")
#     plot(res.mca$ind$coord[,y]+matrix(runif(2*nrow(res.mca$ind$coord[,y]),0,0.1),ncol=2),xlab=xlab,ylab=ylab,col=as.numeric(x@indiv@partition),pch=as.numeric(x@indiv@partition))
#   }
#   
# )