### This algorithm samples the paramters according to their posterior distribution

CMM_bayesparam <- function(x,model,m,chauffe,nbiter){
  ## initialization
  proba <- CMM_MLE(x,model,m,5,tol=10^(-2))$proba
  x <- resume_data(x,model$sigma,m)
  n <- nrow(x);d <- ncol(x);
  alpha <- list();for (k in 1:model$g){alpha[[k]] <- list();for (j in 1:d){alpha[[k]][[j]] <- list()}} 
  
  
  for (iter in 1:chauffe){
    z <- rowSums(t(apply(proba/rowSums(proba),1,cumsum)) < runif(n))+1
    pi <- rdirichlet(1,table(c(1:model$g,z))+3)
    proba <- matrix(pi,n,model$g,byrow=TRUE)
    for (k in unique(z)){
      for (j in 1:d){
        tmp <- sort(table(x[z==k,j]),decreasing=TRUE)+1/2
        tmp_delta <- as.numeric(names(tmp[1:min(model$ell[k,j],length(tmp))]))
        if (length(tmp_delta)<model$ell[k,j]){au <- 1:m[j];au <- au[-tmp_delta];tmp_delta <- c(tmp_delta,sample(au,model$ell[k,j]-length(tmp_delta)))}
        if (model$ell[k,j]>=length(tmp)){tmp <- c(tmp,rep(1/2,model$ell[k,j]+1-length(tmp)))}else{tmp <- c(tmp[1:model$ell[k,j]],sum(tmp[-c(1:model$ell[k,j])]))}
        alpha[[k]][[j]] <- rdirichlet(1,tmp)
        while (any(alpha[[k]][[j]]<(alpha[[k]][[j]][model$ell[k,j]+1]/(m[j]-model$ell[k,j])))){
          alpha[[k]][[j]] <- rdirichlet(1,tmp)
        }
        tmp_alpha <- rep(alpha[[k]][[j]][model$ell[k,j]+1],m[j])
        tmp_alpha[tmp_delta] <- alpha[[k]][[j]][1:model$ell[k,j]]
        proba[,k] <- proba[,k] * tmp_alpha[x[,j]]
      }
    }    
  }
  
  sauv_pi <- matrix(0,nbiter,model$g)
  sauv_alpha <- list();for (k in 1:model$g){sauv_alpha[[k]] <- list();for (j in 1:d){sauv_alpha[[k]][[j]] <- matrix(0,nbiter,m[j])}}
  sauv_z <- matrix(0,n,model$g)
  for (iter in 1:nbiter){
    z <- rowSums(t(apply(proba/rowSums(proba),1,cumsum)) < runif(n))+1
    
    for (k in 1:model$g){sauv_z[which(z==k),k] <- sauv_z[which(z==k),k]+1}
    
    pi <- rdirichlet(1,table(c(1:model$g,z))+1/2)
    sauv_pi[iter,] <- pi
    
    proba <- matrix(pi,n,model$g,byrow=TRUE)
    for (k in unique(z)){
      for (j in 1:d){
        tmp <- sort(table(x[z==k,j]),decreasing=TRUE)+1/2
        tmp_delta <- as.numeric(names(tmp[1:min(model$ell[k,j],length(tmp))]))
        if (length(tmp_delta)<model$ell[k,j]){au <- 1:m[j];au <- au[-tmp_delta];tmp_delta <- c(tmp_delta,sample(au,model$ell[k,j]-length(tmp_delta)))}
        if (model$ell[k,j]>=length(tmp)){tmp <- c(tmp,rep(1/2,model$ell[k,j]+1-length(tmp)))}else{tmp <- c(tmp[1:model$ell[k,j]],sum(tmp[-c(1:model$ell[k,j])]))}
        alpha[[k]][[j]] <- rdirichlet(1,tmp)
        while (any(alpha[[k]][[j]]<(alpha[[k]][[j]][model$ell[k,j]+1]/(m[j]-model$ell[k,j])))){
          alpha[[k]][[j]] <- rdirichlet(1,tmp)
        }
        tmp_alpha <- rep(alpha[[k]][[j]][model$ell[k,j]+1],m[j])
        tmp_alpha[tmp_delta] <- alpha[[k]][[j]][1:model$ell[k,j]]
        proba[,k] <- proba[,k] * tmp_alpha[x[,j]]
        sauv_alpha[[k]][[j]][iter,] <- tmp_alpha
      }
    }    
  }
  
  return(list(sauv_z=sauv_z,sauv_pi=sauv_pi,sauv_alpha=sauv_alpha))
  
}