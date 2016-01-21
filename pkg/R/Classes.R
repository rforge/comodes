###################################################################################
##' Constructor of [\code{\linkS4class{CoModesData}}] class
##'
##'
##' \describe{
##'   \item{data}{data.frame It contains the observations (rows correspond to individual and column to factor)}
##'   \item{modalities}{numeric. Number of levels for each categorical variables.}
##'   \item{levels}{list List of the levels for each categorical variables.}
##' }
##'
##' @examples
##'   getSlots("CoModesData")
##'
##' @name CoModesData-class
##' @rdname CoModesData-class
##' @exportClass CoModesData
##'
setClass(
  Class="CoModesData",
  representation=representation(
    data="data.frame",
    modalities="numeric",
    levels="list"
  ),
  prototype = prototype(
    data=data.frame(),
    modalities=numeric(),
    levels=list()
  )
)

CoModes_data <- function(data, modalities, levels){
  new("CoModesData", data=data, modalities=modalities, levels=levels)
}


###################################################################################
##' Constructor of [\code{\linkS4class{CoModesModel}}] class
##'
##'
##'  
##' \describe{
##'   \item{nbclasses}{numeric. Number of components.}
##'   \item{sigma}{numeric. Give the block where each variable is affiliated.}
##'   \item{modes}{mqtrix. Number of modes (rows correspond to components and columns to blocks)}
##' }
##'
##' @examples
##'   getSlots("CoModesModel")
##'
##' @name CoModesModel-class
##' @rdname CoModesModel-class
##' @exportClass CoModesModel
##'
setClass(
  Class="CoModesModel",
  representation=representation(
    nbclasses="numeric",
    sigma="numeric",
    modes="matrix"
  ),
  prototype = prototype(
    nbclasses=numeric(),
    sigma=numeric(),
    modes=matrix()
  )
)

CoModes_model <- function(nbclasses, sigma, modes){
  new("CoModesModel", nbclasses=nbclasses, sigma=sigma, modes=modes)
}

###################################################################################
##' Constructor of [\code{\linkS4class{CoModesParam}}] class
##'
##'
##' \describe{
##'   \item{proportions}{numeric. Component proportions.}
##'   \item{alpha}{list. Probabilities of the modes.}
##' }
##'
##' @examples
##'   getSlots("CoModesParam")
##'
##' @name CoModesParam-class
##' @rdname CoModesParam-class
##' @exportClass CoModesParam
##'
setClass(
  Class="CoModesParam",
  representation=representation(
    proportions="numeric",
    alpha="list"
  ),
  prototype = prototype(
    proportions=numeric(),
    alpha=list()
  )
)

CoModes_param <- function(proportions, alpha){
  new("CoModesParam", proportions=proportions, alpha=alpha)
}

###################################################################################
##' Constructor of [\code{\linkS4class{CoModesIndiv}}] class
##'
##'
##'  
##' \describe{
##'   \item{proba}{matrix. Posterior probability per component.}
##'   \item{tik}{matrix. Fuzzy partition.}
##'   \item{partition}{numeric. Hard partition.}
##' }
##'
##' @examples
##'   getSlots("CoModesIndiv")
##'
##' @name CoModesIndiv-class
##' @rdname CoModesIndiv-class
##' @exportClass CoModesIndiv
##'
setClass(
  Class="CoModesIndiv",
  representation=representation(
    proba="matrix",
    tik="matrix",
    partition="numeric"
  ),
  prototype = prototype(
    proba=matrix(),
    tik=matrix(),
    partition=numeric()
  )
)


CoModes_indiv <- function(proba, tik, partition){
  new("CoModesIndiv", proba=proba, tik=tik, partition=partition)
}


###################################################################################
##' Constructor of [\code{\linkS4class{CoModesCriteria}}] class
##'
##'
##'  
##' \describe{
##'   \item{loglike}{numeric. Log-likelihood.}
##'   \item{bic}{numeric. BIC criterion.}
##' }
##'
##' @examples
##'   getSlots("CoModesCriteria")
##'
##' @name CoModesCriteria-class
##' @rdname CoModesCriteria-class
##' @exportClass CoModesCriteria
##'
setClass(
  Class="CoModesCriteria",
  representation=representation(
    loglike="numeric",
    bic="numeric"
  ),
  prototype = prototype(
    loglike=numeric(),
    bic=numeric()
  )
)


CoModes_criteria <- function(loglike, bic){
  new("CoModesCriteria", loglike=loglike, bic=bic)
}

###################################################################################
##' Constructor of [\code{\linkS4class{CoModesResults}}] class
##'
##' This S4 class contains the results from the function \link{CoModescluster}.
##'  
##' \describe{
##'   \item{data}{CoModesData. Information about the data.}
##'   \item{model}{CoModesModel. Information about the model.}
##'   \item{param}{CoModesParam. Information about the parameters.}
##'   \item{indiv}{CoModesIndiv. Information about the individuals.}
##'   \item{criteria}{CoModesCriteria. Information about the criteria.}
##' }
##'
##' @examples
##'   getSlots("CoModesResults")
##'
##' @name CoModesResults-class
##' @rdname CoModesResults-class
##' @exportClass CoModesResults
##'
setClass(
  Class="CoModesResults",
  representation=representation(
    data="CoModesData",
    model="CoModesModel",
    param="CoModesParam",
    indiv="CoModesIndiv",
    criteria="CoModesCriteria"
  )
)

CoModes_res <- function(input){
  data <- CoModes_data(input$data, input$modalities, input$levels)
  colnames(input$model$ell) <- paste("Block",1:max(input$model$sigma))
  rownames(input$model$ell) <- paste("Class",1:input$model$g)
  model <- CoModes_model(input$model$g, input$model$sigma, input$model$ell)

  param <- CoModes_param(input$pi, input$alpha)
  indiv <- CoModes_indiv(input$proba, input$tik ,input$partition)
  criteria <- CoModes_criteria(input$loglike, input$bic)
  return(Alpha_organize(new("CoModesResults",data=data, model=model, param=param, indiv=indiv, criteria=criteria)))
}
