setMethod(
  f="barplot",
  signature = c("CoModes_res"),
  definition = function(height,...) {
    
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    # changing marging
    par(mar = rep(2.4,4))
    # decreasing font size
    par(cex = .75)
    
    if (max(height@model@sigma)>3){
      split.screen(c(ceiling(max(height@model@sigma)/2),2))
    }else{
      split.screen(c(max(height@model@sigma),1))
    #  split.screen(c(1,3))
    }
      
    for (j in 1:max(height@model@sigma)){
      screen(j)
      
      
      pro <- matrix(height@param@alpha[[1]][[j]][1:height@model@modes[1,j],1],ncol=1)
      names_var <- list()
      for (h in 1:sum(height@model@sigma==j)){
        names_var[[h]] <- matrix(height@param@alpha[[1]][[j]][1:height@model@modes[1,j],h+1],ncol=1)
      }
      
      if (height@model@nbclasses>1){
        for (k in 2:height@model@nbclasses){
          if (height@model@modes[k,j]<nrow(pro)){
            for (h in 1:sum(height@model@sigma==j)){
              names_var[[h]] <-  cbind(names_var[[h]],c(as.character(height@param@alpha[[k]][[j]][1:height@model@modes[k,j],1+h]),rep(".",nrow(pro)-height@model@modes[k,j]) ))
            }
            pro <- cbind(pro,c(height@param@alpha[[k]][[j]][1:height@model@modes[k,j],1],rep(0,nrow(pro)-height@model@modes[k,j]) ))
          }else if (height@model@modes[k,j]>nrow(pro)){
            for (h in 1:sum(height@model@sigma==j)){
              names_var[[h]] <-  rbind( names_var[[h]],matrix(".",height@model@modes[k,j]-nrow(pro),ncol(pro)))
            }
            pro <- rbind(pro,matrix(0,height@model@modes[k,j]-nrow(pro),ncol(pro)))
            pro <- cbind(pro,height@param@alpha[[k]][[j]][1:height@model@modes[k,j],1])
            for (h in 1:sum(height@model@sigma==j)){
              names_var[[h]] <-  cbind(names_var[[h]],as.character(height@param@alpha[[k]][[j]][1:height@model@modes[k,j],1+h]))
            }
          }else{
            pro <- cbind(pro,height@param@alpha[[k]][[j]][1:height@model@modes[k,j],1])
            for (h in 1:sum(height@model@sigma==j)){
              names_var[[h]] <-  cbind(names_var[[h]],as.character(height@param@alpha[[k]][[j]][1:height@model@modes[k,j],1+h]))
            }
          }    
        }  
      }
      
      pro <- t(pro)
      rownames(pro) <- paste("c",1:height@model@nbclasses)
      mp <- barplot(pro, beside = TRUE, axisnames = TRUE, names.arg=names(t),ylab="P(x^jh=1|z)", ylim=c(0,min(max(pro)*1.5,1.1)))
      
      mtext(1, at = mp, text = paste("C",1:height@model@nbclasses), line = 0, cex = 0.7)
      
      for (h in 1:sum(height@model@sigma==j))
      mtext(3, at = mp, text = as.vector(t(names_var[[h]])), line = -h, cex = 0.7)
      title(paste(paste(names(height@model@sigma)[which(height@model@sigma==j)],collapse="-"),".",sep=""))
    }
    close.screen(all = TRUE)
    # restore plotting parameters
    par(op)
  }
)