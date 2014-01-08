\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,CoModes_res,numeric-method}
\title{Plotting of the individuals of class [\code{\linkS4class{CoModes_res}}]}
\arguments{
  \item{x}{an object of class
           [\code{\linkS4class{CoModes_res}}]}
  
  \item{y}{vector of size indicating the axes of the correspondence analysis where individuals will be drawn.}
  
  \item{...}{further arguments passed to or from other
             methods}
}
\description{
 Plot of the individuals in a multiple correspondence analysis map, 
 the colors and symbols indicate the class membership estimated by CMM.
}
\examples{
# Loading of the data set Acute
data(acute)

# Discretization of the first variable to obtain categorical variables
acute[acute[,1]<37,1] <- 1
acute[acute[,1]>38,1] <- 3
acute[acute[,1]>3,1] <- 2
acute[,1] <- factor(acute[,1],levels=1:3,labels=c("cold","normal","fiever"))

# Model selection and parameter estimation of CMM
res.CoModes <- CoModescluster(acute[,1:6],2)

# Plot of the individuals in a multiple correspondence analysis map, 
# the colors and symbols indicate the class membership estimated by CMM.
plot(res.CoModes,c(2,3))
}
\seealso{
  \code{\link{plot}}
}
