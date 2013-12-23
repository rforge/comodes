## one algorithm EM for estimate the paramters
library(MCMCpack)
library(ade4)

XEM <- function(x,g,ell,m,tol){
  #notations
  n <- nrow(x);  d <- ncol(x) ; proba <- matrix(1/g,n,g)
  # initialization des paramètres
  for (k in 1:g){ for (j in 1:d){alpha_tmp <- rdirichlet(1,rep(1,m[j]));proba[,k] <- proba[,k]*alpha_tmp[x[,j]]}}
  #calcule vraisemblance
  prec <- -Inf; loglike <- sum(log(rowSums(proba)))-10^8;
  
  while((loglike-prec)>tol){
    # E step  
    tik <- sweep(proba,1,rowSums(proba),"/")
    # M step
    proba <- matrix(colSums(tik)/n,n,g,byrow=TRUE);
    for (k in 1:g){
      for (j in 1:d){
        tmp <- rep(0,length(unique(x[,j])))
        names(tmp) <- unique(x[,j])
        loc <- 1
        for (h in unique(x[,j])){
          tmp[loc] <- sum((x[,j]==h)*tik[,k])/sum(tik[,k])
          loc=loc+1
        }
        tmp <- sort(tmp,decreasing=TRUE)
        alpha_tmp <- rep(sum(tmp[-c(1:ell[k,j])])/(m[j]-ell[k,j]),m[j])
        alpha_tmp[as.numeric(names(tmp))[1:ell[k,j]]] <- tmp[1:ell[k,j]]
        proba[,k] <- proba[,k]* alpha_tmp[x[,j]]
      }
    }
    prec <- loglike
    loglike <- sum(log(rowSums(proba)));
    
    # teste fonctionnement EM (à supprimer)
    if (prec>loglike){print("pb dans EM");print(c(prec,loglike))}
  }
  return(list(loglike=loglike,proba=proba))
}





## This function returns the MLE by performing different initialization of EM algorithm
CMM_MLE <- function(x,model,m,nbinit,tol){
  x <- resume_data(x,model$sigma,m)
  m_tmp <- m
  m <- rep(1,max(model$sigma));for (j in 1:max(model$sigma)){m[j] <- prod(m_tmp[which(model$sigma==j)])}
  ref <- XEM(x,model$g,model$ell,m,tol)
  for (it in 2:nbinit){cand <- XEM(x,model$g,model$ell,m,tol); if (cand$loglike>ref$loglike){ref <- cand;}}
  tik <- ref$proba/rowSums(ref$proba)
  pi <- colSums(tik)/sum(tik)
  alpha <- list();
  for (k in 1:model$g){
    alpha[[k]] <- list();
    for (j in 1:ncol(x)){
      tmp <- rep(0,length(unique(x[,j])))
      names(tmp) <- unique(x[,j])
      loc <- 1
      for (h in unique(x[,j])){
        tmp[loc] <- sum((x[,j]==h)*tik[,k])/sum(tik[,k])
        loc=loc+1
      }
      tmp <- sort(tmp,decreasing=TRUE)
      alpha[[k]][[j]] <- rep(sum(tmp[-c(1:model$ell[k,j])])/(m[j]-model$ell[k,j]),model$ell[k,j]+1)
      alpha[[k]][[j]][1:model$ell[k,j]] <- tmp[1:model$ell[k,j]]
      names(alpha[[k]][[j]]) <- c(names(tmp)[1:model$ell[k,j]],"other")
    
    }
  }
  cl <- rep(0,nrow(tik))
  for (k in 1:model$g){cl[which(rowSums(sweep(ref$proba,1,ref$proba[,k],"<="))==model$g)]<-k}
  return(list(model=model,tik=tik,proba=ref$proba,partition=cl,alpha=alpha,pi=pi,loglike=ref$loglike,bic=ref$loglike-(sum(model$ell)+model$g-1)*0.5*log(nrow(x))))
}
