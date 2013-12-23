setClass(
  Class="CoModes_data",
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
  new("CoModes_data", data=data, modalities=modalities, levels=levels)
}



setClass(
  Class="CoModes_model",
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
  new("CoModes_model", nbclasses=nbclasses, sigma=sigma, modes=modes)
}


setClass(
  Class="CoModes_param",
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
  new("CoModes_param", proportions=proportions, alpha=alpha)
}

setClass(
  Class="CoModes_indiv",
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
  new("CoModes_indiv", proba=proba, tik=tik, partition=partition)
}

setClass(
  Class="CoModes_criteria",
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
  new("CoModes_criteria", loglike=loglike, bic=bic)
}


setClass(
  Class="CoModes_res",
  representation=representation(
    data="CoModes_data",
    model="CoModes_model",
    param="CoModes_param",
    indiv="CoModes_indiv",
    criteria="CoModes_criteria"
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
  
  output <- new("CoModes_res",data=data, model=model, param=param, indiv=indiv, criteria=criteria)
  
  return(Alpha_organize(output))
#  return(output)
}
