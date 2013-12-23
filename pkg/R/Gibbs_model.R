### This algorithm sample the models according to their posterior distribution


CMM_model <- function(x,modalites,g,chauffe=500,itermax=2*10^3){
  d <- ncol(x)
  y <- x
  n <- nrow(x)
  m <- modalites
  model <- list(g=g,sigma=1:d,ell=matrix(m-1,g,d,byrow=TRUE))
  proba <- CMM_MLE(x,model,m,1,tol=Inf)$proba
  
  sauv_sigma <- matrix(0,1,d)
  eff_sigma <- 0
  sauv_ell <- list()
  
  for (it in 1:(chauffe+itermax)){

#  print("échantillonnage des appartenances des individus aux classes")
  z <- rowSums(t(apply(proba/rowSums(proba),1,cumsum)) < runif(n))+1
  
  
  ### saut de model
  
  #print("generation du candidat")
    model_cand <- model;var <- sample(1:ncol(x),1);bto <- sample(1:d,1);bfrom <- model$sigma[var];model_cand$sigma[var] <- bto;if (bto==bfrom){model_cand$sigma[var] <- d+1;model_cand$ell<-cbind(model_cand$ell,rep(1,model_cand$g))}
    while(length(unique(model_cand$sigma))<3){model_cand <- model;var <- sample(1:ncol(x),1);bto <- sample(1:d,1);bfrom <- model$sigma[var];model_cand$sigma[var] <- bto;if (bto==bfrom){model_cand$sigma[var] <- d+1;model_cand$ell<-cbind(model_cand$ell,rep(1,model_cand$g))}}  
      mu <- 0;
      
      if (any(model_cand$sigma==bto)){
        m_tmp <- as.numeric(modalites[which(model_cand$sigma==bto)]);  y_tmp <- compute_y(as.matrix(x[,which(model_cand$sigma==bto)]-1),m_tmp)
        mu <- mu - log(length(unique(model_cand$sigma))*(length(unique(model_cand$sigma))>3) + sum(table(model_cand$sigma)>1)*(length(unique(model_cand$sigma))==3)) - model$g*log(prod(m_tmp)-1);
        for (k in 1:model$g){if(any(z==k)){model_cand$ell[k,bto] <- sample(1:(prod(m_tmp)-1),1);    mu <- mu + log_pxz(sort(table(y_tmp[which(z==k)]),decreasing=TRUE),prod(m_tmp),model_cand$ell[k,bto])}}
      }
      
      if (any(model_cand$sigma==bfrom)){
        m_tmp <- as.numeric(modalites[which(model_cand$sigma==bfrom)]); y_tmp <- compute_y(as.matrix(x[,which(model_cand$sigma==bfrom)]-1),m_tmp)

        mu <- mu - log(length(unique(model_cand$sigma))*(length(unique(model_cand$sigma))>3) + sum(table(model_cand$sigma)>1)*(length(unique(model_cand$sigma))==3)   )- model$g*log(prod(m_tmp)-1) ;
        for (k in 1:model$g){model_cand$ell[k,bfrom] <- sample(1:(prod(m_tmp)-1),1);mu <- mu + log_pxz(sort(table(y_tmp[which(z==k)]),decreasing=TRUE),prod(m_tmp),model_cand$ell[k,bfrom])}
      }
      
      for (k in 1:model$g){
        mu <- mu - log_pxz(sort(table(y[which(z==k),bfrom]),decreasing=TRUE),m[bfrom],model$ell[k,bfrom]) + log(length(unique(model$sigma))*(length(unique(model$sigma))>3) + sum(table(model$sigma)>1)*(length(unique(model$sigma))==3)   )+ model$g*log(m[bfrom]-1);
        if (bto!=bfrom){#ai
          mu <- mu -  log_pxz(sort(table(y[which(z==k),bto]),decreasing=TRUE),m[bto],model$ell[k,bto])+  log(length(unique(model$sigma))*(length(unique(model$sigma))>3) + sum(table(model$sigma)>1)*(length(unique(model$sigma))==3)   )+ model$g*log(m[bto]-1);
        }
      } 
      
   #   print("acceptation rejet")
      if (runif(1)<exp(mu)){
        loc <- 0
        d <- length(unique(model_cand$sigma))
        m <- rep(0,d)
        model$ell <- matrix(0,model$g,d)
        for (h in unique(model_cand$sigma)){
          loc <- loc+1
          model$sigma[which(model_cand$sigma==h)] <- loc
          m[loc] <- prod(modalites[which(model_cand$sigma==h)])
          model$ell[,loc] <- model_cand$ell[,h]
        }
        y <- resume_data(x,model$sigma,modalites)
      }
    #}

   proba <- matrix(rdirichlet(1,table(c(1:model$g,z))+3),n,model$g,byrow=TRUE)
  # print(model$ell)
    for (k in 1:model$g){
      for (j in 1:max(model$sigma)){
        if (any(z==k)){
          tmp <- sort(table(y[which(z==k),j]),decreasing=TRUE)
          if (m[j]>2){model$ell[k,j] <- rell_bis(tmp,m[j],model$ell[k,j])}else{model$ell[k,j]=1}
        }else{
          tmp <- 0;
        }
        tmp_delta <- as.numeric(names(tmp[1:min(model$ell[k,j],length(tmp))]))
        tmp <- c(tmp[1:min(length(tmp),model$ell[k,j])] + 0.5 ,sum(tmp) - sum(tmp[1:min(length(tmp),model$ell[k,j])])+0.5)
        if (length(tmp)<(model$ell[k,j]+1)){tmp <- c(tmp,rep(0.5,model$ell[k,j]+1-length(tmp)))}
        al <- rdirichlet(1,tmp)
        cp <- 0
        while (any(al<(al[model$ell[k,j]+1]/(m[j]-model$ell[k,j])))){al <- rdirichlet(1,tmp);cp <- cp+1; if(cp>6){al <- c(tmp[1:model$ell[k,j]],tmp[model$ell[k,j]]/2);al<-al/sum(al)}}
        tmp_alpha <- rep(al[model$ell[k,j]+1],m[j])
        tmp_alpha[tmp_delta] <- al[1:length(tmp_delta)]
        proba[,k] <- proba[,k] * tmp_alpha[y[,j]]
      }
    }
    if (it>chauffe){
      if (any(rowSums(sweep(sauv_sigma,2,model$sigma,"=="))==ncol(x))){
        repere <- which(rowSums(sweep(sauv_sigma,2,model$sigma,"=="))==ncol(x))
        eff_sigma[repere] <- 1 + eff_sigma[repere]
        for (k in 1:model$g){
          for (j in 1:max(model$sigma)){
            if(any(sauv_ell[[repere]][[k]][[j]][,1]==model$ell[k,j])){
              sauv_ell[[repere]][[k]][[j]][which(sauv_ell[[repere]][[k]][[j]][,1]==model$ell[k,j]),2] <- 1 + sauv_ell[[repere]][[k]][[j]][which(sauv_ell[[repere]][[k]][[j]][,1]==model$ell[k,j]),2]
            }else{
              sauv_ell[[repere]][[k]][[j]] <- rbind(sauv_ell[[repere]][[k]][[j]],c(model$ell[k,j],1))
            }
          }
        }
  
      }else{
        eff_sigma <- c(eff_sigma,1)
        repere <- length(eff_sigma)
        sauv_sigma <- rbind(sauv_sigma,model$sigma)
        sauv_ell[[repere]] <- list()
        for (k in 1:model$g){
          sauv_ell[[repere]][[k]] <- list()
          for (j in 1:max(model$sigma)){
            sauv_ell[[repere]][[k]][[j]] <- matrix(c(model$ell[k,j],1),1,2)
          }
        }
      }
      
    }
    
  
  }
#  print(cbind(sauv_sigma,eff_sigma))
  return(list(sauv_sigma=sauv_sigma,eff_sigma=eff_sigma,sauv_modes=sauv_ell)) 
}