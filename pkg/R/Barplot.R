setMethod(
  f="barplot",
  signature = c("CoModes_res"),
  definition = function(height,...) {
    
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
        #  split.screen(c(1,3))
      }
      
      for (j in which(colSums(height@model@modes)>0)){
        screen(j)
        
        if (any(height@model@modes[,j]>0)){
          pro <- height@param@alpha[[1]][[j]][-nrow(height@param@alpha[[1]][[j]]),c(2:ncol(height@param@alpha[[1]][[j]]),1)]
          
          if (height@model@nbclasses>1){
            pro <- cbind(pro,matrix(0,nrow(pro),height@model@nbclasses-1))
            colnames(pro)[-(1:sum(height@model@sigma==j))] <- paste("c",1:height@model@nbclasses,sep="")
            
            for (k in 2:height@model@nbclasses){
              pro[,sum(height@model@sigma==j)+k] <-  height@param@alpha[[k]][[j]][nrow(height@param@alpha[[k]][[j]]),1]
              for (h in 1:(nrow(height@param@alpha[[k]][[j]])-1)){
                if (any(rowSums(sweep(x=as.matrix(pro[,1:sum(height@model@sigma==j)]),MARGIN=2,STATS=as.vector(height@param@alpha[[k]][[j]][h,-1]),FUN="=="))==sum(height@model@sigma==j))){
                  who <- which(rowSums(sweep(x=as.matrix(pro[,1:sum(height@model@sigma==j)]),MARGIN=2,STATS=as.vector(height@param@alpha[[k]][[j]][h,-1]),FUN="=="))==sum(height@model@sigma==j))
                  pro[who,sum(height@model@sigma==j)+k] <- height@param@alpha[[k]][[j]][h,1]
                }else{
                  cp <- data.frame(height@param@alpha[[k]][[j]][h,-1])
                  for (k2 in 1:height@model@nbclasses){
                    cp <- cbind(cp, height@param@alpha[[k2]][[j]][nrow(height@param@alpha[[k2]][[j]]),1])
                  }
                  cp[1,sum(height@model@sigma==j)+k] <-  height@param@alpha[[k]][[j]][h,1]
                  colnames(cp) <- colnames(pro)
                  pro <- rbind(pro,cp)
                  
                }
                
              }
              
              
            }
          }
          
          if (nrow(pro)==1){
            num <- t(as.matrix(pro[,-c(1:sum(height@model@sigma==j))]))
            colnames(num)=NULL
            par(mar=c(sum(height@model@sigma==j),2.5,2.5,2.5))
            mp <- barplot(num)
          }else{
            
            
            num <- t(as.matrix(pro[,-(1:sum(height@model@sigma==j))]))
            or <- order(colSums(num),decreasing=T)
            colnames(num) <- NULL
            par(mar=c(sum(height@model@sigma==j),2.5,2.5,2.5))
            mp <- barplot(num[,or])
          }  
          for (h in 1:sum(height@model@sigma==j))
            mtext(1, at = mp, text = pro[or,h], line = h-1, cex = 0.7)
          
          title(paste(paste(names(height@model@sigma)[which(height@model@sigma==j)],collapse="-"),".",sep=""))
        }
        
        
      }
      close.screen(all = TRUE)
      # restore plotting parameters
      par(op)
    }
  }
)