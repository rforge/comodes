
CMM_ell <- function(x,modalites,model,chauffe=10^3,itermax=3*10^3){
  
  y <- resume_data(x,model$sigma,modalites)
  d <- ncol(y)
  n <- nrow(x)
  m <- rep(0,d);for(h in 1:d){m[h]=prod(modalites[which(model$sigma==h)])}
 # model$ell <- matrix(1,model$g,d,byrow=TRUE)
  proba <- CMM_MLE(x,model,modalites,1,tol=Inf)$proba

  sauv_ell <- list();for (k in 1:model$g){       sauv_ell[[k]] <- list();for( b in 1:d){sauv_ell[[k]][[b]]<- matrix(0,1,2)}}
  
  for (it in 1:(chauffe+itermax)){
    print(it)
    
    ### échantillonnage des appartenances des individus aux classes
    z <- rowSums(t(apply(proba/rowSums(proba),1,cumsum)) < runif(n))+1
    
    
    ### saut de model  
    # echantillonnage des nombres de modes
    
    for (k in 1:model$g){
      for (j in 1:d){
        if (m[j]>2){
          model$ell[k,j] <- rell_bis(sort(table(y[which(z==k),j]),decreasing=TRUE),m[j],model$ell[k,j])
        }else{
          model$ell[k,j] <- 1
        }

      }
    }
    
    ### generation d'un paramètres pour les proba a posteriori
    proba <- matrix(rdirichlet(1,table(c(1:model$g,z))+3),n,model$g,byrow=TRUE)
    for (k in unique(z)){
      for (j in 1:d){
        tmp <- sort(table(y[z==k,j]),decreasing=TRUE)+1/2
        tmp_delta <- as.numeric(names(tmp[1:min(model$ell[k,j],length(tmp))]))
        if (length(tmp_delta)<model$ell[k,j]){au <- 1:m[j];au <- au[-tmp_delta];tmp_delta <- c(tmp_delta,sample(au,model$ell[k,j]-length(tmp_delta)))}
        if (model$ell[k,j]>=length(tmp)){tmp <- c(tmp,rep(1/2,model$ell[k,j]+1-length(tmp)))}else{tmp <- c(tmp[1:model$ell[k,j]],sum(tmp[-c(1:model$ell[k,j])]))}
        al <- rdirichlet(1,tmp)
        cp <- 0
        while (any(al<(al[model$ell[k,j]+1]/(m[j]-model$ell[k,j])))){al <- rdirichlet(1,tmp);cp <- cp+1; if(cp>25){al <- rep(1/m[j],m[j]);print("pb dirichlet")}}
        tmp_alpha <- rep(al[model$ell[k,j]+1],m[j])
        tmp_alpha[tmp_delta] <- al[1:model$ell[k,j]]
        proba[,k] <- proba[,k] * tmp_alpha[y[,j]]
      }
    }
    
    if (it>chauffe){
        for (k in 1:model$g){
          for (j in 1:max(model$sigma)){
            if(any(sauv_ell[[k]][[j]][,1]==model$ell[k,j])){
              sauv_ell[[k]][[j]][which(sauv_ell[[k]][[j]][,1]==model$ell[k,j]),2] <- 1 + sauv_ell[[k]][[j]][which(sauv_ell[[k]][[j]][,1]==model$ell[k,j]),2]
            }else{
              sauv_ell[[k]][[j]] <- rbind(sauv_ell[[k]][[j]],c(model$ell[k,j],1))
            }
          }
        }
    }
    
  }
  ell <- matrix(0,model$g,d)
  for (k in 1:model$g){
    for (j in 1:d){
      if (nrow(sauv_ell[[k]][[j]])>1){
        ell[k,j] <- sauv_ell[[k]][[j]][which(sauv_ell[[k]][[j]][,2]==max(sauv_ell[[k]][[j]][,2]))[1],1]
      }

    }
  }
  return(ell) 
}