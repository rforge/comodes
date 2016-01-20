
CheckInputs <- function(x, g, Gibbs_init,  Gibbs_iter, burnin, EM_init, EM_tol){
  if ( is.data.frame(x) == FALSE)  stop("The data set has to be a data.frame!") 
  if ( any(is.na(x)) == TRUE) stop("Missing values are not allowed!")
  for (j in 1:ncol(x)) {if (class(x[,j]) != "factor")  stop("Each variable must be a factor!")}
  moda <- rep(0,ncol(x))
  levels <- list()
  cp_x <- x
  for (j in 1:ncol(x)){
    levels[[j]] <- levels(x[,j])
    if ( length(levels[[j]])<2)
      stop("All the factors have to be at least two levels!")
    x[,j] <- as.numeric(x[,j])
    moda[j] <- length(levels[[j]])
    
  }
  if ((g != ceiling(g))||(g<1)) stop("Class number has to be integer!")
  if ((Gibbs_iter != ceiling(Gibbs_iter))||(Gibbs_iter<1)) stop("Gibbs_iter has to be integer!")
  if ((burnin != ceiling(burnin))||(burnin<1)) stop("burnin has to be integer!")
  if ((Gibbs_init != ceiling(Gibbs_init))||(Gibbs_init<1)) stop("Gibbs_init has to be integer!")
  if ((EM_init != ceiling(EM_init))||(EM_init<1)) stop("EM_init has to be integer!")
  return(list(moda=moda, levels=levels, g=g, x=x, Gibbs_iter=Gibbs_iter,  burnin=burnin, EM_init=EM_init, EM_tol=EM_tol))
}

ComputeLevelBlock<- function(x,modalites){
  if ((class(x)=="factor") || (class(x)=="numeric")){
    output <- as.numeric(x)
  }else if (class(x)=="data.frame"){
    base <- c(1,exp(cumsum(log(modalites)))[-length(modalites)])
    for(j in 1:length(modalites)) x[,j] <- as.numeric(as.character(x[,j]))
    x <- as.matrix(x, ncol=length(modalites))
    output <- (x-1)%*%base+1
  }else{
    stop("error in the function ComputeLevelBlock")
  }
  return(output)
}

ComputeLevelAllBlock <- function(x, sigma, modalites){
  y <- matrix(0,nrow(x),max(sigma))
  for (j in 1:max(sigma)){
    if (sum(sigma==j)>1)
      y[,j] <- ComputeLevelBlock(x[,which(sigma==j)],modalites[which(sigma==j)])
    else
      y[,j] <- x[,which(sigma==j)]
  }
  return(y)
}



LogIntComDataLikeBlockCompo <- function(ech, m, el){
  px <- 0
  if (length(ech)<=el){ech <- c(ech,rep(0,el-length(ech)+1))}
  px <- -sum(ech[(el+1):length(ech)])*log(m-el);
  if (el>0){
    for (tmp in 1:el)    px <- px +Ibeta(1/(m-tmp+1),ech[tmp]+1,sum(ech[(tmp+1):length(ech)])+1,lower=FALSE,log=TRUE) - log(m-tmp) 
  }  
  return(px)
}


CondSamplingModeNumber <- function(ech, currentmodel, k, b){
  currentmodes <- currentmodel$ell[k,b]
  ratio <- 1
  if (currentmodel$ell[k,b] == 0){
    candidatesmodes <- currentmodel$ell[k,b] + 1
  }else if (currentmodel$ell[k,b] == (currentmodel$blocklevels[b]-1)){
    candidatesmodes <- currentmodel$ell[k,b]- 1
  }else{
    candidatesmodes <- currentmodel$ell[k,b] + sample(c(-1,1),1)
  }
  if ( (candidatesmodes==1) || (candidatesmodes==(currentmodel$blocklevels[b]-1))) ratio <- ratio / 2
  if ( (currentmodes==1) || (currentmodes==(currentmodel$blocklevels[b]-1))) ratio <- ratio * 2
  ratio <- log(ratio) + LogIntComDataLikeBlockCompo(ech, currentmodel$blocklevels[b], candidatesmodes) - LogIntComDataLikeBlockCompo(ech, currentmodel$blocklevels[b], currentmodes)
  if (runif(1)<exp(ratio)) currentmodel$ell[k,b] <- candidatesmodes
  return(currentmodel)
}

xfromy <- function(y,m){
  if (length(m)>1){
    y <- y-1
    for (h in 2:length(m)) m[h] <- m[h] * m[h-1]
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
          for (h in 2:ncol(output@param@alpha[[k]][[j]])){
            output@param@alpha[[k]][[j]][,h] <- as.character(c(output@data@levels[[vbles[h-1]]][output@param@alpha[[k]][[j]][1:output@model@modes[k,j],h]],"."))
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
