
XEM_latent <- function(x,g,tol=10^(-3),m){
  n <- nrow(x);d <- ncol(x);
  proba <- matrix(1/g,n,g);pi <- rep(1/g,g);
  for (k in 1:g){
    for (j in 1:d){
      tmp_alpha <- rdirichlet(1,rep(1,m[j]))
      proba[,k] <- proba[,k]*tmp_alpha[x[,j]]
    }    
  }
  prec <- -Inf;
  loglike <- sum(log(rowSums(proba)))
  
  while((loglike-prec)>tol){
    tik <- sweep(proba,1,rowSums(proba),"/")
    
    pi <- colSums(tik)/n
    proba <- matrix(pi,n,g,byrow=TRUE)
    for (k in 1:g){
      for (j in 1:d){
        tmp_alpha <- rep(0,m[j])
        for (h in 1:m[j]){ tmp_alpha[h] <- sum((x[,j]==h)*tik[,k])/sum(tik[,k])}
        proba[,k] <- proba[,k]*tmp_alpha[x[,j]]
      }
    }
    prec <- loglike;
    loglike <- sum(log(rowSums(proba)))
  }
  return(proba)
}


CIMclustering <- function(x,g,nbinit=50,tol=10^(-3)){
  m <- rep(0,ncol(x));for (h in 1:ncol(x)){m[h] <- length(unique(x[,h]))}
 ref=list(loglike=-Inf)
 for (it in 1:nbinit){
   cand <- XEM_latent(x,g,tol,m)
   if (sum(log(rowSums(cand)))>ref$loglike){
     ref <- list(proba=cand,loglike=sum(log(rowSums(cand))))
   }
 }
 
 ref$nbparam <- (g-1) + g*sum((m-1))
 ref$tik <- sweep(ref$proba,1,rowSums(ref$proba),"/")
 ref$alpha <- list();
 for (j in 1:ncol(x)){
   ref$alpha[[j]] <- matrix(0,g,m[j])
   for (k in 1:g){
     tmp_alpha <- rep(0,m[j])
     for (h in 1:m[j]){ tmp_alpha[h] <- sum((x[,j]==h)*ref$tik[,k])/sum(ref$tik[,k])}
     ref$alpha[[j]][k,] <- tmp_alpha
   }
 }
 ref$pi <- colSums(ref$tik)/nrow(x)
 ref$cl <- rep(0,nrow(x))
 for (k in 1:g){ref$cl[which(rowSums(sweep(ref$proba,1,ref$proba[,k],"<="))==g)]<-k}
 ref$bic=ref$loglike - ref$nbparam*0.5*log(nrow(x))
  return(ref)
}

