

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


# compute_y <- function(x,modalites){
#   if (length(modalites)>1){
#     modalites <- c(1,modalites[-length(modalites)])
#     for (loc in 2:length(modalites)){modalites[loc] <- modalites[loc]*modalites[loc-1]} 
#   }
#   return((x%*%modalites)+1)
# }
# 
# resume_data <- function(x,sigma,modalites){
#   y <- matrix(0,nrow(x),max(sigma))
#   for (j in 1:max(sigma)){
#     if (sum(sigma==j)>1){
#       y[,j] <- compute_y(as.matrix(x[,which(sigma==j)])-1,modalites[which(sigma==j)])
#     }else{
#       y[,j] <- x[,which(sigma==j)]
#     }
#   }
#   return(y)
# }
# 
# log_pxz <- function(ech,m,el){
#   px <- 0
#   if (length(ech)<=el){ech <- c(ech,rep(0,el-length(ech)+1))}
#   px <- -sum(ech[(el+1):length(ech)])*log(m-el);
#   if (el>0){
#     for (tmp in 1:el)    px <- px +Ibeta(1/(m-tmp+1),ech[tmp]+1,sum(ech[(tmp+1):length(ech)])+1,lower=FALSE,log=TRUE) - log(m-tmp) 
#   }  
#   return(px)
# }
# 
# rell_bis<- function(ech,m,actu){
#   if ((actu>1)&&(actu<(m-1))){
#     cand <- actu-1+2*(runif(1)<1/2);
#     if (any(c(1,m-1)==cand)){coeff <- 2}else{coeff <- 1}
#     cand <- sort(c(actu,cand));
#   }else if (actu==1){cand <- c(1,2);coeff <- 0.5;}else{cand <- c(actu-1,actu);coeff <- 0.5;}
#   if (length(ech)<=cand[2]){ech <- c(ech,rep(0,cand[2]))}
#     px <- c(-sum(ech[(cand[1]+1):length(ech)])*log(m-cand[1]) ,-sum(ech[(cand[2]+1):length(ech)])*log(m-cand[2]) )
#     px[2] <- px[2] -log(m-cand[2]) +  Ibeta(1/(m-cand[2]+1),ech[cand[2]]+1,sum(ech[(cand[2]+1):length(ech)])+1,lower=FALSE,log=TRUE)
#     px <- exp(px - max(px))
#     px <- coeff*px/sum(px)
#   return(sample(cand,1,prob=px))
# }