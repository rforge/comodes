\docType{methods}
\name{barplot}
\alias{barplot}
\alias{barplot,CoModes_res-method}
\title{Plotting of the parameters of class [\code{\linkS4class{CoModes_res}}]}
\arguments{
  \item{height}{an object of class [\code{\linkS4class{CoModes_res}}]}
  
  
  \item{...}{further arguments passed to or from other
             methods}
}
\description{
Plot the parameters of each component
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

# Plot of the parameters
barplot(res.CoModes)}
\seealso{
  \code{\link{barplot}}
}
