\name{CoModescluster}
\alias{CoModescluster}
\title{Create an instance of the [\code{\linkS4class{CoModes_res}}] class}
\usage{
  CoModescluster(x, g, Gibbs_init=2,  Gibbs_iter=50,
  Gibbs_chauffe=50, EM_init=5, EM_tol=10^(-3) )
}
\arguments{
  \item{x}{frame containing qualitative data. Rows correspond to individuals and
  columns correspond to variables. Each columns must be a factor.}

  \item{g}{the number of classes.}

  \item{Gibbs_init}{the number of different initializations of  a Gibbs chain
  performing the model selection.}

  \item{Gibbs_iter}{the number of iterations of a Gibbs chain
  performing the model selection.}

  \item{Gibbs_chauffe}{the number of first iterations of a Gibbs chain
  which are not taken into account in the estimation of the posterior distribution of the 
  model.}

  \item{EM_init}{the number of different initializations of  the EM algorithm
  performing the maximum likelihood estimation for the model selected by the Gibbs chain.}

  \item{EM_tol}{stopping criteria. The EM algorithm is stopped when the difference between to successive log likelihood is smaller than EM_tol.}
}
\value{
  Returns an instance of the
  [\code{\linkS4class{CoModes_res}}] class composed by five attributes: \describe{
  \item{data}{an instance of the [\code{\linkS4class{CoModes_data}}]  class containing the elements relative to the data.}
  \item{model}{an instance of the [\code{\linkS4class{CoModes_model}}]  class containing the elements relative to the model.}
  \item{param}{an instance of the [\code{\linkS4class{CoModes_param}}]  class containing the elements relative to the parameters.}
  \item{indiv}{an instance of the [\code{\linkS4class{CoModes_indiv}}]  class containing the elements relative to the individuals.}
  \item{criteria}{an instance of the [\code{\linkS4class{CoModes_criteria}}]  class containing the elements relative to the criteria.}
    
    }
}
\description{
  This function computes an optimal mixture model according
  to the integrated likelihood and computes its associated
  maximum likelihood estimates.
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
