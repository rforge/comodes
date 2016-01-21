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
##' # Data Loading
##' data(alzheimer)
##' # Model selection and Parameter estimation for CMM with 2 components (8 MCMC chains are used for model selection)
##' results <- CoModescluster(alzheimer, 2, 8)
##' # Summary of the results
##' summary(results)
##' # Display the probabilities of the modes
##' barplot(results)
##' }
##' 

NULL

##' Real categorical data: acute
##' 
##' The data was created by a medical expert as a data set to test the expert system, which will perform the presumptive diagnosis of two diseases of urinary system. The basis for rules detection was Rough Sets Theory. Each instance represents an potential patient. Here the first variable has been modified to only take integer value.
##'
##' Informations are available here: \url{http://archive.ics.uci.edu/ml/datasets/Acute+Inflammations}
##' 
##' \itemize{
##' \item{obs: }{Patients described by six variables Temperature of patient (35C-42C), occurrence of nausea (yes,no),
##'  lumbar pain (yes,no), urine pushing (continuous need for urination) (yes,no), mMicturition pains (yes,no), burning of urethra, itch, swelling of urethra outlet (yes,no)}
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
##' \item{obs: }{Observed variables}
##' \item{class: }{normal find, metastases, malign lymph, fibrosis}
##' }
##'
##' Informations are available here: \url{http://archive.ics.uci.edu/ml/datasets/Lymphography}
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
##' \item{obs: }{plumage and external morphological characteristics}
##' \item{class: }{dichrous, lherminieri and subalaris}
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

##' Real categorical data: alzheimer
##' 
##' Presence or absence of 6 symptoms of Alzheimer's disease (AD) in 240 patients diagnosed with early onset AD conducted in the Mercer Institute in St. James's Hospital, Dublin.
##' A binary matrix, consisting of 240 rows and 6 columns, with each row denoting an individual and each column denoting the presence/absence of one of the 6 symptoms: Hallucination, Activity, Aggression, Agitation, Diurnal and Affective. A 1 denotes the presence of a symptom, a 0 the absence
##' 
##' @references Moran M, Walsh C, Lynch A, Coen RF, Coakley D and Lawlor B (2004). Syndromes of behavioural and psychological symptoms in mild Alzheimer's disease. International Journal of Geriatric Psychiatry, 19(4), 359-364.
##' @references White A and Muprhy T (2014). BayesLCA: An R Pacakge for Bayesian Latent Class Analysis. Journal of Statistical Software, 61(13),1-28                                   
##' @name alzheimer
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(alzheimer)
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


###################################################################################
##' Create an instance of the [\code{\linkS4class{CoModesResults}}] class
##'
##' This function performs the model selection and the parameter inference.
##' x, g, Gibbs_init=30,  Gibbs_iter=min(4000,(ncol(x)*400)), burnin=min(ncol(x)*400,4000), EM_init=25, EM_tol=10^(-3), nbcores=Gibbs_init
##' @param x data.frame, where each column is a factor.
##' @param g integer, defines the number of components.
##' @param Gibbs_init integer, number of chains performed for model selection (default 30).
##' @param Gibbs_iter integer, number of iterations of the MCMC algorithm for model selection (default min(ncol(x)*400,4000))
##' @param burnin integer, number of iterations of the burn-in of the MCMC algorithm (default min(ncol(x)*400,4000))
##' @param EM_init integer, number of runs of EM algorithm for parameter inference (default 25)
##' @param EM_tol numeric, tolerance for the stopping criterion of the EM algorithm (default 0.001)
##' @param nbcores number of cores used by the algorithm (only for Linux and MAC). (default Gibbs_init)
##'
##' @examples
##' \dontrun{
##' # Data Loading
##' data(alzheimer)
##' # Model selection and Parameter estimation for CMM with 2 components (8 MCMC chains are used for model selection)
##' results <- CoModescluster(alzheimer, 2, 8)
##' # Summary of the results
##' summary(results)
##' # Display the probabilities of the modes
##' barplot(results)
##' }
##'
##' @return Returns an instance of the [\code{\linkS4class{CoModesResults}}] class. 
##' @export
##'
##'
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
