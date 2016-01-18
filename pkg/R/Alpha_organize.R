xfromy <- function(y,m){
  if (length(m)>1){
    y <- y-1
    m <- as.integer(exp(cumsum(log(m))))
    x <- matrix(0,length(y),length(m))
    for (h in length(m):2){
      x[,h] <- y %/% m[h-1]
      y <- y - (y%/%m[h-1])*m[h-1]
    }
    x[,1] <- y
  }else{
    x <- y
  }
  
  return(x+1)
}

Alpha_organize <- function(output){
  for (k in 1:output@model@nbclasses){
    for (j in unique(output@model@sigma)){
      vbles <- which(output@model@sigma==j)
      m <- output@data@modalities[vbles]
      a <- output@param@alpha[[k]][[j]]
      if (output@model@modes[k,j]>0){
        if (length(vbles)>1){
          output@param@alpha[[k]][[j]] <- cbind(as.numeric(output@param@alpha[[k]][[j]]),rbind(xfromy(as.numeric(names(a[-length(a)])),m),rep(NA,length(m))))        
          output@param@alpha[[k]][[j]] <- data.frame(output@param@alpha[[k]][[j]])
          who <- which(output@model@sigma == j)
          for (h in 2:ncol(output@param@alpha[[k]][[j]])){
            output@param@alpha[[k]][[j]][,h] <- as.character(c(output@data@levels[[who[h-1]]][as.numeric(output@param@alpha[[k]][[j]][1:(nrow(output@param@alpha[[k]][[j]])-1),h])],"."))
          }
          output@param@alpha[[k]][[j]][,1] <- as.numeric(output@param@alpha[[k]][[j]][,1])
        }else{
          output@param@alpha[[k]][[j]] <-   data.frame(probability=as.numeric(a),tmp=as.character(c(output@data@levels[[vbles]][as.numeric(names(a)[-length(a)])],".")))
        }
      }else{
        output@param@alpha[[k]][[j]] <- data.frame(as.numeric(a))
        output@param@alpha[[k]][[j]] <- cbind(output@param@alpha[[k]][[j]], matrix(".", 1, length(vbles)))
      }
            colnames(output@param@alpha[[k]][[j]]) <- c("probability",names(output@model@sigma[which(output@model@sigma==j)]))
    }
  }
  return(output)
}