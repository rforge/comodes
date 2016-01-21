###################################################################################
##' Barplot of a class [\code{\linkS4class{CoModesResults}}]  
##'
##' Barplot of qualitative data from a [\code{\linkS4class{CoModesResults}}] object using parameters
##' to plot probablities of modalities.
##'
##' Each line corresponds to one block of variable. Barplot is drawn for each cluster with the probabilities for 
##' each level (with is a mode for at least on cluster) to be in that cluster.
##'  
##' @param height an object of class [\code{\linkS4class{CoModesResults}}]
#' @name barplot
#' @rdname barplot-methods
#' @docType methods
#' @exportMethod summary
#' 
#' 
NULL

#' @rdname barplot-methods
#' @aliases barplot barplot,CoModesResults-method
#'

setMethod(
  f="barplot",
  signature = c("CoModesResults"),
  definition = function(height) {
    
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    # changing marging
    par(mar = rep(2.4,4))
    # decreasing font size
    par(cex = .75)
    
    nbv <- sum(colSums(height@model@modes)>0)
    if (nbv>0){
      if (nbv>3){
        split.screen(c(ceiling(nbv/2),2))
      }else{
        split.screen(c(nbv,1))
      }
      
      for (j in which(colSums(height@model@modes)>0)){
        screen(j)
        vari <- which(height@model@sigma==j)
        if (length(vari) == 1){
          modes <- levels(height@param@alpha[[1]][[j]][,-1])
          for (k in 2:height@model@nbclasses) modes <- c(modes, levels(height@param@alpha[[k]][[j]][,-1]))
          modes <- unique(modes)
          proba <- matrix(NA, length(modes),height@model@nbclasses)
          for (k in 1:height@model@nbclasses){
            for(h in 1:length(modes)){
              if (any(levels(height@param@alpha[[k]][[j]][,2])==modes[h])){
                proba[h,k] <- height@param@alpha[[k]][[j]][which(height@param@alpha[[k]][[j]][,2]==modes[h]),1]
              }else{
                proba[h,k] <- height@param@alpha[[k]][[j]][height@model@modes[k,j]+1,1]
              }
            }
          }
          if (any(modes==".")) modes[which(modes==".")] <- "Other"
        }else{
          modes <- height@param@alpha[[1]][[j]][,-1]
          for (k in 2:height@model@nbclasses) modes <- rbind(modes, height@param@alpha[[k]][[j]][,-1])
          modes <- as.matrix(modes)
          for (h in 1:nrow(modes)){
            if (h<=nrow(modes))
            who <- which(rowSums(sweep(modes, 2, modes[h,], "==")) == ncol(modes))
            if (length(who)>1){
              modes <- modes[- who[which(who!=h)],]
            } 
          }
          proba <- matrix(NA, nrow(modes),height@model@nbclasses)
          for (k in 1:height@model@nbclasses){
            for(h in 1:nrow(modes)){
              if (any(apply(sweep(height@param@alpha[[k]][[j]][,-1], 2, modes[h,], "=="),1,"all"))){
                proba[h,k] <- height@param@alpha[[k]][[j]][which(apply(sweep(height@param@alpha[[k]][[j]][,-1], 2, modes[h,], "=="),1,"all")),1]
              }else{
                proba[h,k] <- height@param@alpha[[k]][[j]][height@model@modes[k,j]+1,1]
              }
            }
          }
          if (any(modes[,1]==".")) modes[which(modes[,1]=="."),] <- c("Other",rep("", ncol(modes)-1))
       }
        ord <- order(rowSums(proba), decreasing = TRUE)

       
       
          par(mar=c(sum(height@model@sigma==j),2.5,2.5,2.5))
          mp <-  barplot(t(proba[ord,]))
       
       if (length(vari) == 1){
          for (h in 1:sum(height@model@sigma==j))   mtext(1, at = mp, text = modes[ord], line = h-1, cex = 0.7)
       }else{
         for (h in 1:sum(height@model@sigma==j))   mtext(1, at = mp, text = modes[ord,h], line = h-1, cex = 0.7)
         
       }
          title(paste(paste(names(height@model@sigma)[which(height@model@sigma==j)],collapse="-"),".",sep=""))
        
        
        
      }
      close.screen(all.screens =  TRUE)
      # restore plotting parameters
      par(op)
    }
  }
)