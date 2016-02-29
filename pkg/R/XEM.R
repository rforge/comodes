## one algorithm EM for estimate the paramters
XEM <- function(x, model, tol){
  #notations
  n <- nrow(x);  d <- ncol(x) ; proba <- matrix(1/model$g,n,model$g)
  # initialization des paramètres
  for (k in 1:model$g){ for (j in 1:d){alpha_tmp <- rdirichlet(1,rep(1,model$blocklevels[j]));proba[,k] <- proba[,k]*alpha_tmp[x[,j]]}}
  #calcule vraisemblance
  prec <- -Inf; loglike <- sum(log(rowSums(proba)))-10^8;  
  cp <-1
  while((loglike-prec)>tol){
    cp <- cp + 1
    # E step  
    tik <- sweep(proba,1,rowSums(proba),"/")
    # M step
    proba <- matrix(log(colSums(tik)/sum(tik)), n, model$g, byrow=TRUE);
    for (k in 1:model$g){
      for (j in 1:d){
        if (model$ell[k,j]>0){
          tmp <- rep(0, model$blocklevels[j])
          names(tmp) <- 1:length(tmp)
          for (h in unique(x[,j]))  tmp[h] <- sum(tik[which(x[,j]==h),k])
          tmp <- sort(tmp,decreasing=TRUE)
          modes <- as.numeric(names(tmp))[1:model$ell[k,j]]
          tmp <- c(tmp[1:model$ell[k,j]]+1, sum(tmp)-sum(tmp[1:model$ell[k,j]])+1) / (sum(tik[,k]) +model$ell[k,j]+1)
          alpha_tmp <- rep(tmp[model$ell[k,j]+1], model$blocklevels[j])/(model$blocklevels[j]-model$ell[k,j])
          alpha_tmp[modes] <- tmp[1:model$ell[k,j]]
        }else{
          alpha_tmp <- rep(1/model$blocklevels[j],model$blocklevels[j])
        }
        proba[,k] <- proba[,k] + log(alpha_tmp[x[,j]])
      }
    }
    maxlogproba <- apply(proba,1,max)
    proba <- exp(sweep(proba, 1, maxlogproba,"-"))
    prec <- loglike
    loglike <- sum(maxlogproba) + sum(log(rowSums(proba)));
  }
  tik=sweep(proba,1,rowSums(proba),"/")
  alpha <- list()
  for (k in 1:model$g){
    alpha[[k]]<-list()
    for (j in 1:d){
      if (model$ell[k,j]>0){
        tmp <- rep(0, model$blocklevels[j])
        names(tmp) <- 1:length(tmp)
        for (h in unique(x[,j]))  tmp[h] <- sum(tik[which(x[,j]==h),k])
        tmp <- sort(tmp,decreasing=TRUE)
        modes <- as.numeric(names(tmp))[1:model$ell[k,j]]
        tmp <- c(tmp[1:model$ell[k,j]]+1, sum(tmp)-sum(tmp[1:model$ell[k,j]])+1) / (sum(tik[,k])+model$ell[k,j]+1)
        alpha[[k]][[j]] <- tmp
        alpha[[k]][[j]][model$ell[k,j]+1] <- alpha[[k]][[j]][model$ell[k,j]+1]/(model$blocklevels[j] - model$ell[k,j])
        names(alpha[[k]][[j]])  <- c(modes,"other")
      }else{
        alpha[[k]][[j]] <- rep(1/model$blocklevels[j],1)
        names(alpha[[k]][[j]]) <- "other"
      }
    }
  }
  return(list(loglike=loglike,proba=proba,  tik=tik, pi= colSums(tik)/sum(tik), alpha=alpha, cl=apply(tik,1,which.max) ))
}

## This function returns the MLE by performing different initialization of EM algorithm
CMM_MLE <- function(x, model, nbinit, tol){
  ref <- XEM(x,model,tol)
  for (it in 2:nbinit){cand <- XEM(x,model,tol); if (cand$loglike>ref$loglike){ref <- cand;}}
  return(list(model=model, tik=ref$tik, proba=ref$proba, partition=ref$cl, alpha=ref$alpha, pi=ref$pi, loglike=ref$loglike,bic=ref$loglike-(sum(model$ell)+model$g-1)*0.5*log(nrow(x))))
}

