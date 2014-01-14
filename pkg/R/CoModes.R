#'
#' 
NULL 



#' CoModescluster function
#' 
#' This function performs clustering for categorical data using the block model extension of the latent class model.
#'  
#' 
#' 
#' 
#' 
#' @export
#' 
#' @name CoModescluster
#' @useDynLib CoModes
#'

CoModescluster <- function(x, g, Gibbs_init=2,  Gibbs_iter=50, Gibbs_chauffe=50, EM_init=5, EM_tol=10^(-3) ){
  if ( is.data.frame(x) == FALSE){
    stop("The data set has to be a data.frame!")
  }
  
  if ( any(is.na(x)) == TRUE){
    stop("Missing values are not allowed!")
  }
  
  if ((g != ceiling(g))||(g<1))
    stop("Class number has to be integer!")
    
  
  moda <- rep(0,ncol(x))
  levels <- list()
  cp_x <- x
  for (j in 1:ncol(x)){
    levels[[j]] <- levels(x[,j])
    if ( length(levels[[j]])<2)
      stop("All the factors have to be at least two levels!")
      
    x[,j] <- as.numeric(x[,j])
    moda[j] <- length(unique(x[,j]))
  }
  

  output <- list(bic=-Inf)
  for (ch in 1:Gibbs_init){
    res <- CMM_model(x, moda, g, Gibbs_chauffe, Gibbs_iter)
    best_sigma <- res$sauv_sigma[which(res$eff_sigma==max(res$eff_sigma))[1],]
    tmp_ell <- res$sauv_modes[[which(res$eff_sigma==max(res$eff_sigma))[1]]]
    best_ell <- matrix(0,g,max(best_sigma))
    for (k in 1:g){
      for (j in 1:max(best_sigma)){
        best_ell[k,j] <- tmp_ell[[k]][[j]][which(tmp_ell[[k]][[j]][,2]==max(tmp_ell[[k]][[j]][,2]))[1],1]
      }
    }
    
    res_model <- CMM_MLE(x, list(g=g, sigma=best_sigma, ell=best_ell), moda, EM_init, EM_tol)
    if ( res_model$bic > output$bic )
      output <- res_model
    
  }
  
 
  output$data <- cp_x
  output$modalities <- moda
  output$levels <- levels
  names(output$model$sigma) <- colnames(x)
 # print(output$alpha)
  output <- CoModes_res(output)
  
  return(output)
}

