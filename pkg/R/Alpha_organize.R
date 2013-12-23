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
      m <- output@data@modalities[which(output@model@sigma == j)]
      a <- output@param@alpha[[k]][[j]]
      if (sum(output@model@sigma==j)>1){
        output@param@alpha[[k]][[j]] <- cbind(as.numeric(output@param@alpha[[k]][[j]]),rbind(xfromy(as.numeric(names(a[-length(a)])),m),rep(NA,length(m))))        
        output@param@alpha[[k]][[j]] <- data.frame(output@param@alpha[[k]][[j]])
        who <- which(output@model@sigma == j)
        for (h in 2:ncol(output@param@alpha[[k]][[j]])){
          output@param@alpha[[k]][[j]][,h] <- as.character(c(output@data@levels[[who[h-1]]][as.numeric(output@param@alpha[[k]][[j]][1:(nrow(output@param@alpha[[k]][[j]])-1),h])],"."))
        }
        output@param@alpha[[k]][[j]][,1] <- as.numeric(output@param@alpha[[k]][[j]][,1])
       
      }else{
        output@param@alpha[[k]][[j]] <- data.frame(probability=as.numeric(output@param@alpha[[k]][[j]]),tmp=as.character(c(output@data@levels[[which(output@model@sigma == j)]][as.numeric(names(a)[-length(a)])],".")))

      }
      colnames(output@param@alpha[[k]][[j]]) <- c("probability",names(output@model@sigma[which(output@model@sigma==j)]))

    }
  }
  return(output)
}