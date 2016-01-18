## one algorithm EM for estimate the paramters
XEM <- function(x, model, tol){
  #notations
  n <- nrow(x);  d <- ncol(x) ; proba <- matrix(1/model$g,n,model$g)
  # initialization des paramètres
  for (k in 1:model$g){ for (j in 1:d){alpha_tmp <- rdirichlet(1,rep(1,model$blocklevels[j]));proba[,k] <- proba[,k]*alpha_tmp[x[,j]]}}
  #calcule vraisemblance
  prec <- -Inf; loglike <- sum(log(rowSums(proba)))-10^8;  
  while((loglike-prec)>tol){
    # E step  
    tik <- sweep(proba,1,rowSums(proba),"/")
    # M step
    proba <- matrix(colSums(tik)/n,n,model$g,byrow=TRUE);
    for (k in 1:model$g){
      for (j in 1:d){
        if (model$ell[k,j]>0){
          tmp <- rep(0, model$blocklevels[j])
          names(tmp) <- 1:length(tmp)
          for (h in unique(x[,j]))  tmp[h] <- sum((x[,j]==h)*tik[,k])/sum(tik[,k])
          tmp <- sort(tmp,decreasing=TRUE)
          alpha_tmp <- rep(sum(tmp[-c(1:model$ell[k,j])])/(model$blocklevels[j]-model$ell[k,j]),model$blocklevels[j])
          alpha_tmp[as.numeric(names(tmp))[1:model$ell[k,j]]] <- tmp[1:model$ell[k,j]]
        }else{
          alpha_tmp <- rep(1/model$blocklevels[j],model$blocklevels[j])
        }
        proba[,k] <- proba[,k]* alpha_tmp[x[,j]]
      }
    }
    prec <- loglike
    loglike <- sum(log(rowSums(proba)));
  }
  return(list(loglike=loglike,proba=proba))
}
## This function returns the MLE by performing different initialization of EM algorithm
CMM_MLE <- function(x, model, nbinit, tol){
  ref <- XEM(x,model,tol)
  for (it in 2:nbinit){cand <- XEM(x,model,tol); if (cand$loglike>ref$loglike){ref <- cand;}}
  tik <- ref$proba/rowSums(ref$proba)
  pi <- colSums(tik)/sum(tik)
  alpha <- list();
  for (k in 1:model$g){
    alpha[[k]] <- list();
    for (j in 1:ncol(x)){
      if (model$ell[k,j]>0){
        tmp <- rep(0, model$blocklevels[j])
        names(tmp) <- 1:length(tmp)
        for (h in unique(x[,j])) tmp[h] <- sum((x[,j]==h)*tik[,k])/sum(tik[,k])
        tmp <- sort(tmp,decreasing=TRUE)
        alpha[[k]][[j]] <- rep(sum(tmp[-c(1:model$ell[k,j])])/(model$blocklevels[j]-model$ell[k,j]),model$ell[k,j]+1)
        alpha[[k]][[j]][1:model$ell[k,j]] <- tmp[1:model$ell[k,j]]
        names(alpha[[k]][[j]]) <- c(names(tmp)[1:model$ell[k,j]],"other")
      }else{
        alpha[[k]][[j]] <- 1/model$blocklevels[j]
        names(alpha[[k]][[j]]) <- "other"
      }
      
    }
  }
  cl <- apply(tik, 1, which.max)
  return(list(model=model,tik=tik,proba=ref$proba,partition=cl,alpha=alpha,pi=pi,loglike=ref$loglike,bic=ref$loglike-(sum(model$ell)+model$g-1)*0.5*log(nrow(x))))
}

