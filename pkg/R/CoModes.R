#'
#' 
NULL 

CheckInputs <- function(x, g, Gibbs_init,  Gibbs_iter, burnin, EM_init, EM_tol){
  if ( is.data.frame(x) == FALSE)  stop("The data set has to be a data.frame!") 
  if ( any(is.na(x)) == TRUE) stop("Missing values are not allowed!")
  for (j in 1:ncol(x)) {if (class(x[,j]) != "factor")  stop("Each variable must be a factor!")}
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
  if ((g != ceiling(g))||(g<1)) stop("Class number has to be integer!")
  if ((Gibbs_iter != ceiling(Gibbs_iter))||(Gibbs_iter<1)) stop("Gibbs_iter has to be integer!")
  if ((burnin != ceiling(burnin))||(burnin<1)) stop("burnin has to be integer!")
  if ((Gibbs_init != ceiling(Gibbs_init))||(Gibbs_init<1)) stop("Gibbs_init has to be integer!")
  if ((EM_init != ceiling(EM_init))||(EM_init<1)) stop("EM_init has to be integer!")
  return(list(moda=moda, levels=levels, g=g, x=x, Gibbs_iter=Gibbs_iter,  burnin=burnin, EM_init=EM_init, EM_tol=EM_tol))
}

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

OneAnalysis <- function(nb, input){
  res <- CMM_model(input$x, input$moda, input$g, input$burnin, input$Gibbs_iter)
  selectedmodel <- list(g=input$g, sigma=res$sauv_sigma[which(res$eff_sigma==max(res$eff_sigma))[1],])
  tmp_ell <- res$sauv_modes[[which(res$eff_sigma==max(res$eff_sigma))[1]]]
  selectedmodel$ell <- matrix(0, input$g, max(selectedmodel$sigma))
  for (k in 1:input$g){
    for (j in 1:max(selectedmodel$sigma))  selectedmodel$ell[k,j] <- tmp_ell[[k]][[j]][which(tmp_ell[[k]][[j]][,2]==max(tmp_ell[[k]][[j]][,2]))[1],1]
  }
  selectedmodel$blocklevels <- rep(0, max(selectedmodel$sigma))
  for (j in 1:max(selectedmodel$sigma)) selectedmodel$blocklevels[j] <- prod(input$moda[which(selectedmodel$sigma==j)])
  y <- ComputeLevelAllBlock(input$x, selectedmodel$sigma, input$moda)
  res_model <- CMM_MLE(y, selectedmodel,  input$EM_init, input$EM_tol)
  return(res_model)
}

CoModescluster <- function(x, g, Gibbs_init=5,  Gibbs_iter=10**3, burnin=10**3, EM_init=25, EM_tol=10^(-3), nbcores=Gibbs_init){
  input <- CheckInputs(x, g, Gibbs_init,  Gibbs_iter, burnin, EM_init, EM_tol)
  nbcores <- min(detectCores(all.tests = FALSE, logical = FALSE), nbcores)
  if ((nbcores>1)&&(Sys.info()["sysname"] != "Windows")){
    output <- mclapply(X = as.list(1:Gibbs_init),
                       FUN = OneAnalysis,
                       input=input,
                       mc.cores = nbcores, mc.preschedule = TRUE, mc.cleanup = TRUE)
  }else{
    output <- list()
    for (it in 1:Gibbs_init) output[[it]] <- OneAnalysis(it, input)
  }  
  allbic <- rep(-Inf, Gibbs_init)
  for (it in 1:Gibbs_init) allbic <- output[[it]]$bic
  output <- output[[which.max(allbic)]]
  
  output$data <- x
  output$modalities <- input$moda
  output$levels <- input$levels
  names(output$model$sigma) <- colnames(x)
  output <- CoModes_res(output)
  return(output)
}