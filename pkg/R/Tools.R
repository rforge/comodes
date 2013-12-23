### functions tools
library(zipfR)


### This function converts the modalities crossings in a new categorical variable
compute_y <- function(x,modalites){
  if (length(modalites)>1){
    modalites <- c(1,modalites[-length(modalites)])
    for (loc in 2:length(modalites)){modalites[loc] <- modalites[loc]*modalites[loc-1]} 
  }
  return((x%*%modalites)+1)
}

resume_data <- function(x,sigma,modalites){
  y <- matrix(0,nrow(x),max(sigma))
  for (j in 1:max(sigma)){
    if (sum(sigma==j)>1){
      y[,j] <- compute_y(as.matrix(x[,which(sigma==j)])-1,modalites[which(sigma==j)])
    }else{
      y[,j] <- x[,which(sigma==j)]
    }
  }
  return(y)
}


### Computate of the integrate complete-data log-likelihood
# log_pxz<- function(ech,m,el){
#   px <- 0
#   if (length(ech)<=el){ech <- c(ech,rep(0,el-length(ech)+1))}
#   px <- -sum(ech[(el+1):length(ech)])*log(m-el);
#   for (tmp in 1:el){
#     px <- px + Ibeta(1/(m-tmp+1),ech[tmp]+1/2,sum(ech[(tmp+1):length(ech)]+1/2),lower=FALSE,log=TRUE) - log(m-tmp+1) - Ibeta(1/(m-tmp+1),1/2,1/2,lower=FALSE,log=TRUE)
#   }
#   return(px)
# }

log_pxz<- function(ech,m,el){
  px <- 0
  if (length(ech)<=el){ech <- c(ech,rep(0,el-length(ech)+1))}
  px <- -sum(ech[(el+1):length(ech)])*log(m-el);
  for (tmp in 1:el){
    px <- px +Ibeta(1/(m-tmp+1),ech[tmp]+1,sum(ech[(tmp+1):length(ech)])+1,lower=FALSE,log=TRUE) - log(m-tmp) 
  }
  return(px)
}

# 
# nb_modes_exact<- function(ech,m){
#   n <- sum(ech)
#   px <- rep(-Inf,m)
#   px[1] <- -n*log(m)
#   for (el in 1:(length(ech)-1)){
#     px[el+1] <- -sum(ech[(el+1):length(ech)])*log(m-el);
#     for (tmp in 1:el){
#       px[el+1] <- px[el+1] + Ibeta(1/(m-tmp+1),ech[tmp]+1,sum(ech[(tmp+1):length(ech)]+1),lower=FALSE,log=TRUE) - log(m-tmp) }
#     if ((el>4)&& (all(px[el+1]<px[(el-2):el]))){
#       break;
#     }
#   }
#   return(which(px==max(px))[1]-1)
# }

rell_bis<- function(ech,m,actu){
  if ((actu>1)&&(actu<(m-1))){
    cand <- actu-1+2*(runif(1)<1/2);
    if (any(c(1,m-1)==cand)){coeff <- 2}else{coeff <- 1}
    cand <- sort(c(actu,cand));
  }else if (actu==1){cand <- c(1,2);coeff <- 0.5;}else{cand <- c(actu-1,actu);coeff <- 0.5;}
  if (length(ech)<=cand[2]){ech <- c(ech,rep(0,cand[2]))}
    px <- c(-sum(ech[(cand[1]+1):length(ech)])*log(m-cand[1]) ,-sum(ech[(cand[2]+1):length(ech)])*log(m-cand[2]) )
    px[2] <- px[2] -log(m-cand[2]) +  Ibeta(1/(m-cand[2]+1),ech[cand[2]]+1,sum(ech[(cand[2]+1):length(ech)])+1,lower=FALSE,log=TRUE)
    px <- exp(px - max(px))
    px <- coeff*px/sum(px)
  return(sample(cand,1,prob=px))
}

# rell<- function(ech,m){
#   if (length(ech)<m){ech <- c(ech,rep(0,m-length(ech)))}
#   n <- sum(ech)
#   px <- rep(-n*log(m),m)
#   for (el in 1:(length(ech)-1)){
#     px[el+1] <- -sum(ech[(el+1):length(ech)])*log(m-el);
#     for (tmp in 1:el){
#       px[el+1] <- px[el+1] + Ibeta(1/(m-tmp+1),ech[tmp]+1/2,sum(ech[(tmp+1):length(ech)]+1/2),lower=FALSE,log=TRUE) - log(m-tmp+1) - Ibeta(1/(m-tmp+1),1/2,1/2,lower=FALSE,log=TRUE)}
#   }
#   px <- px[-1]
#   px <- exp(px - max(px))
#   px <- px/sum(px)
#   return(sample(1:length(px),1,prob=px))
# }
# 
# 
# rell_bis<- function(ech,m,actu){
#   if ((actu>1)&&(actu<(m-1))){
#     cand <- actu-1+2*(runif(1)<1/2);
#     if (any(c(1,m-1)==cand)){coeff <- 2}else{coeff <- 1}
#     cand <- sort(c(actu,cand));
#   }else if (actu==1){cand <- c(1,2);coeff <- 0.5;}else{cand <- c(actu-1,actu);coeff <- 0.5;}
#   if (length(ech)<=cand[2]){ech <- c(ech,rep(0,cand[2]))}
#     px <- c(-sum(ech[(cand[1]+1):length(ech)])*log(m-cand[1]) ,-sum(ech[(cand[2]+1):length(ech)])*log(m-cand[2]) )
#     px[2] <- px[2] -log(m-cand[2]+1) +  Ibeta(1/(m-cand[2]+1),ech[cand[2]]+1/2,sum(ech[(cand[2]+1):length(ech)]+1/2),lower=FALSE,log=TRUE)  - Ibeta(1/(m-cand[2]+1),1/2,1/2,lower=FALSE,log=TRUE)
#     px <- exp(px - max(px))
#     px <- coeff*px/sum(px)
#   return(sample(cand,1,prob=px))
# }
# 




# build_summary <- function(obj,m){
#   moda <- rep(1,max(obj$model$sigma))
#   for (h in 1:length(moda)){moda[h] <- prod(m[which(obj$model$sigma==h)])}
#   kappa <- obj$model$ell
#   for (k in 1:nrow(kappa)){kappa[k,] <- kappa[k,]/(moda-1)}
#   tau <- matrix(0,nrow(kappa),ncol(kappa))
#   for (k in 1:nrow(tau)){
#     for (j in 1:ncol(tau)){
#       tau[k,j] <- sum(obj$alpha[[k]][[j]][1:obj$model$ell[k,j]])
#     }
#   }
#   return(list(kappa=kappa,tau=tau))
# }
# 
# present_alpha <- function(obj,m){
#   resume <- list()
#   for (j in 1:max(obj$model$sigma)){
#     resume[[j]] <- list()
#     for (k in 1:obj$model$g){
#       resume[[j]][[k]] <- cbind(xfromy(as.numeric(names(obj$alpha[[k]][[j]])[1:obj$model$ell[k,j]]),m[which(obj$model$sigma==j)]),obj$alpha[[k]][[j]][1:obj$model$ell[k,j]])
#     }
#   }
#   return(resume)  
# }