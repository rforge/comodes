\name{CoModes_data-class}
\alias{CoModes_data-class}
\title{The [\code{\linkS4class{CoModes_data}}] class}
\description{
  This is a class to contain results from CoModescluster function.
}
\details{
  This class is composed by three attributs
  \describe{
    \item{data}{a data.frame where the rows correspond to the individuals and where the columns correspond to the variables.}
    \item{modalities}{a vector indicating the modality number of each variable.}
    \item{levels}{the elements of this list contain the different modalities of each variable.}
  } 
}
\examples{
  getSlots("CoModes_data")
}
