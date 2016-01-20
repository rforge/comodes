##' CoModes a package for clustering categorical data
##'
##' CoModes is a tool for clustering categorical data. The clustering goal is achieved by a mixture model which considers the
##' intra-class dependencies since the observed variables are grouped into conditionally independent blocks. The block distribution
##' is parsimonious since each block follows a multinomial distribution per modes. Under this distribution, the free parameters
##' correspond to the probabilities of the most probable levels (the modes) while the other levels (the non modes) are 
##' assumed to be equiprobable.
##'
##' \tabular{ll}{
##'   Package: \tab CoModes\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 1.0.0\cr
##'   Date: \tab 2016-01-20\cr 
##'   License: \tab GPL-2\cr 
##'   LazyLoad: \tab yes\cr
##' }
##'
##'
##' @name CoModes-package
##' @aliases CoModes
##' @rdname CoModes-package
##' @docType package
##' @keywords package
##' @import methods
##' @import zipfR
##' @import parallel
##' @import MCMCpack
##' @export CoModescluster
##' @exportMethod summary
##' @exportMethod barplot
##' @exportClass CoModesResults
##'
##' @author
##' Author: Marbac Matthieu, Biernacki Christophe and Vandewalle Vincent
##'
##' @references Marbac Matthieu, Biernacki Christophe and Vandewalle Vincent (2015). Latent class model with conditional dependency per modes to cluster categorical data. arXiv:1402.5103
##'
##' @examples
##' \dontrun{
##' a=1
##' }
##' 

NULL

##' Real categorical data: acute
##' 
##' The data was created by a medical expert as a data set to test the expert system, which will perform the presumptive diagnosis of two diseases of urinary system. The basis for rules detection was Rough Sets Theory. Each instance represents an potential patient. Here the first variable has been modified to only take integer value.
##'
##' Informations are available here: http://archive.ics.uci.edu/ml/datasets/Acute+Inflammations
##' 
##' \itemize{
##' \item{acute[,1]: }{Temperature of patient (35C-42C)}
##' \item{acute[,2]: }{Occurrence of nausea (yes,no)}
##' \item{acute[,3]: }{Lumbar pain (yes,no)}
##' \item{acute[,4]: }{Urine pushing (continuous need for urination) (yes,no)}
##' \item{acute[,5]: }{Micturition pains (yes,no)}
##' \item{acute[,6]: }{Burning of urethra, itch, swelling of urethra outlet (yes,no)}
##' \item{urin: }{Inflammation of urinary bladder (yes,no)}
##' \item{renal: }{Nephritis of renal pelvis origin (yes,no)}
##' \item{urinrenal: }{Inflammation of urinary bladder  and Nephritis of renal pelvis origin (1:no-no, 2:yes-no, 3:no-yes, 4:yes-yes)} 
##'}
##'
##' 
##' @references J.Czerniak, H.Zarzycki, Application of rough sets in the presumptive diagnosis of urinary system diseases, Artifical Inteligence and Security in Computing Systems, ACS'2002 9th International Conference Proceedings, Kluwer Academic Publishers,2003, pp. 41-51 
##' 
##' @name acute
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(acute)
NULL



##' Real categorical data: lymphography
##' 
##' This lymphography domain was obtained from the University Medical Centre, Institute of Oncology, Ljubljana, Yugoslavia.  Thanks go to M. Zwitter and M. Soklic for providing the data.
##'
##' \itemize{
##' \item{x:}{Observed variables}
##' \item{class: }{normal find, metastases, malign lymph, fibrosis}
##' }
##'
##' Informations are available here: http://archive.ics.uci.edu/ml/datasets/Lymphography
##' 
##' 
##' @name lymphography
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(lymphography)
NULL


##' Real categorical data: seabirds
##' 
##' The Seabirds data set is a biological data set describing 153 puffins (seabirds) by five plumage and external morphological characteristics. These seabirds are divided into three subspecies dichrous (84 birds), lherminieri (34 birds) and subalaris (35 birds).
##' 
##' 
##' \itemize{
##' \item{birds:}{plumage and external morphological characteristics}
##' \item{species: }{dichrous, lherminieri and subalaris}
##' }
##' 
##' 
##' @references Bretagnolle, V. (2007). Personal communication. source: Museum.
##' 
##' @name seabirds
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(seabirds)
NULL


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



CoModescluster <- function(x, g, Gibbs_init=30,  Gibbs_iter=min(4000,(ncol(x)*400)), burnin=min(ncol(x)*400,4000), EM_init=25, EM_tol=10^(-3), nbcores=Gibbs_init){
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
