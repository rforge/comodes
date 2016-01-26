GenerateCandidate <- function(model, modalities){
  model_cand <- model
  # sampling of the blocks
  var <- sample(1:length(model$sigma), 1)
  bto <- sample(1:max(model$sigma), 1)
  bfrom <- model$sigma[var]
  model_cand$sigma[var] <- bto
  if (bto==bfrom){
    bto <- max(model$sigma) + 1
    model_cand$sigma[var] <- max(model$sigma)+1
    model_cand$ell <- cbind(model_cand$ell, rep(1,model_cand$g))
  }
  ratioproba <- max(model$sigma)/max(model_cand$sigma)
  # sampling of the modes  
  model_cand$blocklevels[bto] <- prod(modalities[which(model_cand$sigma==bto)])
  model_cand$ell[,bto] <- sample(1:model_cand$blocklevels[bto], nrow(model_cand$ell), replace=TRUE) - 1
  if (any(model$sigma == bto)){
    ratioproba <- ratioproba * model$blocklevels[bto] / model_cand$blocklevels[bto]
  }else{
    ratioproba <- ratioproba / model_cand$blocklevels[bto]
  }
  if (any(model_cand$sigma == bfrom)){
    model_cand$blocklevels[bfrom] <- prod(modalities[which(model_cand$sigma==bfrom)])
    model_cand$ell[,bfrom] <- sample(1:model_cand$blocklevels[bfrom], nrow(model_cand$ell), replace=TRUE) - 1
    ratioproba <- ratioproba * model$blocklevels[bfrom] / model_cand$blocklevels[bfrom]
  }else{
    model_cand$blocklevels[bfrom] <- 0
    model_cand$ell[,bfrom] <- 0
    ratioproba <- ratioproba * model$blocklevels[bfrom]
  }  
  return(list(model=model_cand, bto=bto, bfrom=bfrom, var=var, ratioproba=ratioproba, modelseed=model))
}

CopyModel <- function(input, modalities){
  output <- input
  output$blocklevels <- rep(0, length(unique(input$sigma)))
  output$ell <- matrix(0, output$g, length(unique(input$sigma)))
  loc <- 0
  for (h in unique(input$sigma)){
    loc <- loc+1
    output$sigma[which(input$sigma==h)] <- loc
    output$blocklevels[loc] <- prod(modalities[which(input$sigma==h)])
    output$ell[,loc] <- input$ell[,h]
  }
  return(output)
}


### This algorithm sample the models according to their posterior distribution
CMM_model <- function(x, modalities, g, burnin, nbiter){
  y <- x
  for (j in 1:ncol(x)) y[,j] <- as.numeric(x[,j])
  currentmodel <- list(g=g, sigma=1:ncol(x), ell=matrix(modalities-1,g, ncol(x),byrow=TRUE), blocklevels=modalities)
  proba <- CMM_MLE(y, currentmodel, 1, 0.1)$proba
  sauv_sigma <- matrix(0, 1, ncol(x))
  eff_sigma <- 0
  sauv_ell <- list()  
  for (it in 1:(burnin+nbiter)){
    #échantillonnage des appartenances des individus aux classes
    z <- apply(proba, 1, function(sa) sample(1:length(sa),1, prob = sa))
    ### saut de model
    # generation du candidat
    candidate <- GenerateCandidate(currentmodel, modalities)
    # ratio contains is the logarithm of the ratio between the probability of the model generation
    ratio <- log(candidate$ratioproba)
    # ratio is updated on the part of the integrated complete-data likelihood of the candidate and current models
    ycandidatebto <- ComputeLevelBlock(x[,which(candidate$model$sigma==candidate$bto)], modalities[which(candidate$model$sigma==candidate$bto)])
    if (any(candidate$model$sigma==candidate$bfrom)) ycandidatebfrom <- ComputeLevelBlock(x[,which(candidate$model$sigma==candidate$bfrom)], modalities[which(candidate$model$sigma==candidate$bfrom)])
    for (k in 1:g){
      if (any(z==k)){
        ratio <- ratio + LogIntComDataLikeBlockCompo(sort(table(ycandidatebto[which(z==k)]),decreasing=TRUE), candidate$model$blocklevels[candidate$bto], candidate$model$ell[k,candidate$bto]) - 
          LogIntComDataLikeBlockCompo(sort(table(y[which(z==k),candidate$bfrom]),decreasing=TRUE), currentmodel$blocklevels[candidate$bfrom], currentmodel$ell[k,candidate$bfrom]) 
        if (any(candidate$model$sigma==candidate$bfrom)) ratio <- ratio + LogIntComDataLikeBlockCompo(sort(table(ycandidatebfrom[which(z==k)]),decreasing=TRUE), candidate$model$blocklevels[candidate$bfrom], candidate$model$ell[k,candidate$bfrom])
        if (any(currentmodel$sigma==candidate$bto)) ratio <- ratio - LogIntComDataLikeBlockCompo(sort(table(y[which(z==k), candidate$bto]), decreasing=TRUE), currentmodel$blocklevels[candidate$bto], currentmodel$ell[k,candidate$bto])
      }      
    }
    # acceptation rejet
    if (runif(1)<exp(ratio)){
      currentmodel <- CopyModel(candidate$model, modalities)
      y <- ComputeLevelAllBlock(x, currentmodel$sigma, modalities)
    }
    

    
    # computation of the conditional probabilities of the component memberships by sampling the continuous parameters of the current model
    proba <- matrix(log(rdirichlet(1,table(c(1:g,z))-1/2)), length(z), g, byrow=TRUE)
    for (k in 1:g){
      if (any(z==k)){
        for (b in 1:max(currentmodel$sigma)){
          # sampling of the number of modes
          tmp <- sort(table(y[which(z==k), b]), decreasing=TRUE)
          currentmodel <- CondSamplingModeNumber(tmp, currentmodel, k, b)
          if (currentmodel$ell[k,b]>0){
            tmp <- sort(table(y[which(z==k), b]), decreasing=TRUE)
            if (length(tmp)< currentmodel$ell[k,b]){
              missing <- currentmodel$ell[k,b] - length(tmp)
              nam <- as.numeric(names(tmp))
              other <- 1:currentmodel$blocklevels[b]
              tmp <- c(tmp, rep(0, missing))
              names(tmp) <- c(nam, sample(other[-nam], missing))
            } 
            eff <- c(tmp[1:currentmodel$ell[k,b]], sum(tmp) - sum(tmp[1:currentmodel$ell[k,b]]))
            tmpdelta <- as.numeric(names(eff[1:currentmodel$ell[k,b]]))
            tmpalpha <- rep(0, length(eff))
            limit <- 0
            while ((limit<6) && (min(tmpalpha[1:currentmodel$ell[k,b]]) <= (tmpalpha[1+currentmodel$ell[k,b]]/(currentmodel$blocklevels[b]-currentmodel$ell[k,b])))){
              tmpalpha <- rdirichlet(1,eff+1)
              limit <- limit + 1
            }
            if (limit==6) tmpalpha <- rep(1/currentmodel$blocklevels[b],currentmodel$blocklevels[b])
            tmpalpha[currentmodel$ell[k,b]+1] <- tmpalpha[currentmodel$ell[k,b]+1]/(currentmodel$blocklevels[b]-currentmodel$ell[k,b])
          }else{
            tmpalpha <- 1/currentmodel$blocklevels[b]
          }
          tmpproba <- rep(tmpalpha[1+currentmodel$ell[k,b]], currentmodel$blocklevels[b])
          if (currentmodel$ell[k,b]>0) tmpproba[tmpdelta] <- tmpalpha[1:currentmodel$ell[k,b]]
          proba[,k] <- proba[,k] + log(tmpproba[y[,b]])
        }
      }
    } 
    proba <- exp(sweep(proba,1,apply(proba,1,max),"-"))
    if (it>burnin){
      repere <- which(rowSums(sweep(sauv_sigma,2,currentmodel$sigma,"=="))==ncol(x))
      if (length(repere)>0){        
        eff_sigma[repere] <- 1 + eff_sigma[repere]
        for (k in 1:currentmodel$g){
          for (j in 1:max(currentmodel$sigma)) sauv_ell[[repere]][[j]][currentmodel$ell[k,j]+1 ,k] <- sauv_ell[[repere]][[j]][currentmodel$ell[k,j]+1 ,k]  + 1
        }        
      }else{
        eff_sigma <- c(eff_sigma,1)
        sauv_sigma <- rbind(sauv_sigma,currentmodel$sigma)
        repere <- length(eff_sigma)
        sauv_ell[[repere]] <- list()
        for (j in 1:max(currentmodel$sigma)) {
          sauv_ell[[repere]][[j]] <- matrix(0, currentmodel$blocklevels[j]+1, currentmodel$g)
          for (k in 1:currentmodel$g) sauv_ell[[repere]][[j]][currentmodel$ell[k,j]+1 ,k] <- 1
        }
      }
    }   
  }
  return(list(sauv_sigma=sauv_sigma,eff_sigma=eff_sigma,sauv_modes=sauv_ell)) 
}
