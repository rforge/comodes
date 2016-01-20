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
    ratioproba <- ratioproba * (model$blocklevels[bto]-1) / (model_cand$blocklevels[bto]-1)
  }else{
    ratioproba <- ratioproba / (model_cand$blocklevels[bto]-1)
  }
  if (any(model_cand$sigma == bfrom)){
    model_cand$blocklevels[bfrom] <- prod(modalities[which(model_cand$sigma==bfrom)])
    model_cand$ell[,bfrom] <- sample(1:model_cand$blocklevels[bfrom], nrow(model_cand$ell), replace=TRUE) - 1
    ratioproba <- ratioproba * (model$blocklevels[bfrom]-1) / (model_cand$blocklevels[bfrom]-1)
  }else{
    model_cand$blocklevels[bfrom] <- 0
    model_cand$ell[,bfrom] <- 0
    ratioproba <- ratioproba * (model$blocklevels[bfrom]-1)
  }  
  return(list(model=model_cand, bto=bto, bfrom=bfrom, var=var, ratioproba=ratioproba))
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
    #échantillonnage des appartenances des individus aux classes"
    z <- apply(proba, 1, function(sa) sample(1:length(sa),1, prob = sa))
    #save(z, currentmodel, file = "current.rda")
    ### saut de model
    # generation du candidat
    candidate <- GenerateCandidate(currentmodel, modalities)
    #save(candidate, file = "candidate.rda")
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
    # sampling of the number of modes
    for (k in 1:g){
      if (any(z==k)){
        for (b in 1:max(currentmodel$sigma)) currentmodel <- CondSamplingModeNumber(sort(table(y[which(z==k), b]), decreasing=TRUE), currentmodel, k, b)
      }
    }
    # computation of the conditional probabilities of the component memberships by sampling the continuous parameters of the current model
    proba <- matrix(rdirichlet(1,table(c(1:g,z))+3), length(z), g, byrow=TRUE)
    for (k in 1:g){
      if (any(z==k)){
        for (b in 1:max(currentmodel$sigma)){
          if (currentmodel$ell[k,b]>0){
            eff <- c(sort(table(c(1:currentmodel$blocklevels[b],y[which(z==k), b])), decreasing=TRUE))
            eff <- c(eff[1:currentmodel$ell[k,b]], sum(eff) - sum(eff[1:currentmodel$ell[k,b]]))
            tmpdelta <- as.numeric(names(eff[1:currentmodel$ell[k,b]]))
            tmpalpha <- rep(0, length(eff))
            limit <- 0
            while ((limit<6) && (min(tmpalpha[1:currentmodel$ell[k,b]]) <= (tmpalpha[1+currentmodel$ell[k,b]]/(currentmodel$blocklevels[b]-currentmodel$ell[k,b])))){
              tmpalpha <- rdirichlet(1,eff)
              limit <- limit + 1
            }
            if (limit==6) tmpalpha <- eff/sum(eff)
            tmpalpha[currentmodel$ell[k,b]+1] <- tmpalpha[currentmodel$ell[k,b]+1]/(currentmodel$blocklevels[b]-currentmodel$ell[k,b])
          }else{
            tmpalpha <- 1/currentmodel$blocklevels[b]
          }
          tmpproba <- rep(tmpalpha[1+currentmodel$ell[k,b]], currentmodel$blocklevels[b])
          if (currentmodel$ell[k,b]>0) tmpproba[tmpdelta] <- tmpalpha[1:currentmodel$ell[k,b]]
          proba[,k] <- proba[,k] * tmpproba[y[,b]]
        }
      }
    }    
    if (it>burnin){
      if (any(rowSums(sweep(sauv_sigma,2,currentmodel$sigma,"=="))==ncol(x))){
        repere <- which(rowSums(sweep(sauv_sigma,2,currentmodel$sigma,"=="))==ncol(x))
        eff_sigma[repere] <- 1 + eff_sigma[repere]
        for (k in 1:currentmodel$g){
          for (j in 1:max(currentmodel$sigma)){
            if(any(sauv_ell[[repere]][[k]][[j]][,1]==currentmodel$ell[k,j])){
              sauv_ell[[repere]][[k]][[j]][which(sauv_ell[[repere]][[k]][[j]][,1]==currentmodel$ell[k,j]),2] <- 1 + sauv_ell[[repere]][[k]][[j]][which(sauv_ell[[repere]][[k]][[j]][,1]==currentmodel$ell[k,j]),2]
            }else{
              sauv_ell[[repere]][[k]][[j]] <- rbind(sauv_ell[[repere]][[k]][[j]],c(currentmodel$ell[k,j],1))
            }
          }
        }        
      }else{
        eff_sigma <- c(eff_sigma,1)
        repere <- length(eff_sigma)
        sauv_sigma <- rbind(sauv_sigma,currentmodel$sigma)
        sauv_ell[[repere]] <- list()
        for (k in 1:currentmodel$g){
          sauv_ell[[repere]][[k]] <- list()
          for (j in 1:max(currentmodel$sigma)) sauv_ell[[repere]][[k]][[j]] <- matrix(c(currentmodel$ell[k,j],1),1,2)
        }
      }
    }   
  }
  return(list(sauv_sigma=sauv_sigma,eff_sigma=eff_sigma,sauv_modes=sauv_ell)) 
}
