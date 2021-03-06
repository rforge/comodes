###################################################################################
##' Summary of a class [\code{\linkS4class{CoModesResults}}]  
##'
##'  
##' @param object an object of class [\code{\linkS4class{CoModesResults}}]
##'
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' 
#' 
NULL

#' @rdname summary-methods
#' @aliases summary summary,CoModesResults-method
#'

setMethod(
  f="summary",
  signature = c("CoModesResults"),
  definition = function(object) {

    cat("**************************************************************************************\n\n")
    cat("Number of variables:", ncol(object@data@data), "    Number of individuals:",length(object@indiv@partition), "\n" )
    
    
    cat("Number of modalities:",object@data@modalities,  "\n" )
    
    
    cat("Class number:", object@model@nbclasses , "    log-likelihood:",object@criteria@loglike,  "    BIC:",object@criteria@bic, "\n")
    
    cat("Block of the variables:", object@model@sigma, "\n")
    
    cat("\n*************************************\n")
    
    
    
    moda <- rep(1,max(object@model@sigma))
    nam <- rep(NA,length(moda))
    for (h in 1:length(moda)){
      moda[h] <- prod(object@data@modalities[which(object@model@sigma==h)])
      nam[h] <- paste(names(object@model@sigma)[which(object@model@sigma==h)],collapse="-")
    }
    kappa <- object@model@modes
    colnames(kappa) <- nam
    
    cat("Mode number:","\n")
    print(kappa)
    
    cat("\n\n**************************************************************************************\n\n")
    for (k in 1:nrow(kappa)){kappa[k,] <- kappa[k,]/(moda-1)}
    tau <- matrix(0,nrow(kappa),ncol(kappa))
    for (k in 1:nrow(tau)){
      for (j in 1:ncol(tau)){
        if (object@model@modes[k,j]>0) tau[k,j] <- sum(object@param@alpha[[k]][[j]][1:object@model@modes[k,j],1])
      }
    }
    
    rownames(tau) <- rownames(kappa)
    colnames(tau) <- colnames(kappa)
    
    cat("AlphakjDot index:", "\n" )
    print(round(tau,2))
    
    cat("\n*************************************\n\n")
    
    cat("Ukjbar index", "\n")
    print(round(kappa,2))
    cat("\n**************************************************************************************\n")
  }
)