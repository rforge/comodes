\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,CoModes_res-method}
\title{Produce result summaries of a CoModes_res class}
\arguments{
  \item{x}{Either a \code{\linkS4class{CoModes_res}} object.}
  
  \item{...}{further arguments passed to or from other
             methods}
}
\value{
  NULL. Summaries to standard out.
}
\description{
  Produce result summaries of a CoModes_res class
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

# Summary of the model
summary(res.CoModes)
}
\seealso{
  \code{\link{summary}}
}
