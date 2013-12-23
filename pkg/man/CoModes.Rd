\docType{package}
\name{CoModes-package}
\alias{CoModes}
\alias{CoModes-package}
\title{Conditional MODES mixture models}
\description{
  CoModes is a tool for fitting, on a categorical data set,
  a mixture model of conditional modes distributions,
  which are multinational distribution where few free parameters correspond to the
  largest probability locations.
}
\details{
  \tabular{ll}{ Package: \tab CoModes\cr Type: \tab
  Package\cr Version: \tab 1.0.0\cr Date: \tab
  2013-12-19\cr License: \tab GPL \cr
  LazyLoad: \tab yes\cr }

  CoModes focuses on the clustering for categorical data. It also provides a partition 
  of the individuals and parameters to characterize each class. Thus, CoModes gives a 
  summary of a complex data set by extracting the main information from the data.
  
  CoModes solves clustering problems by estimating densities with specific multinomial
  distribution, named Conditional Multinomial Model. Under this model, variables are 
  grouped into conditional independent blocks. The corresponding block distribution is a parsimonious multinomial distribution where the few free parameters correspond to the most likely modality crossings, while the remaining probability mass is uniformly spread over the other modality crossings. Thus, the proposed model allows to bring out the intra-class dependence between variables and to summarize each class by a few characteristic modality crossings.
  
  The exact computation of the integrated complete-data log-likelihood allows to perform the model selection, by a Gibbs sampler, reducing the computing time consuming by parameter estimation and avoiding BIC criterion biases pointed out by our experiments. Finally, the parameter estimation is only performed for the best model by an EM algorithm.
  
  The function CoModescluster allows to perform the model selection and  the maximum likelihood estimation. It returns an object of CoModes_res. CoModes contains three tools functions barplot, plot and summary, respectively giving a graphical parameter summary,
  a scatter-plot of the individuals in a multiple correspondence analysis map with a color and a symbol indicating a the class membership, and a summary of the parameters.
  
}



\examples{
# Package loading
require(CoModes)

# Loading of the data set Acute
data(acute)

# Discretization of the first variable to obtain categorical variables
acute[acute[,1]<37,1] <- 1
acute[acute[,1]>38,1] <- 3
acute[acute[,1]>3,1] <- 2
acute[,1] <- factor(acute[,1],levels=1:3,labels=c("cold","normal","fiever"))

# Model selection and parameter estimation of CMM
res.CoModes <- CoModescluster(acute[,1:6],2)

# Plot of the parameters
barplot(res.CoModes)

# Plot of the individuals in a multiple correspondance analysis map, 
# the colors and symbols indicate the class membership estimated by CMM.
plot(res.CoModes,c(2,3))

# Summary of the model
summary(res.CoModes)

}

\author{
  Matthieu Marbac and Christophe Biernacki and Vincent Vandewalle.
}
\references{
  Marbac M., Biernacki C., Vandewalle V., 2014. 
  "Finite mixture model of conditional dependencies modes to cluster categorical data".  Reprint.
}
