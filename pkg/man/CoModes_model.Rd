\name{CoModes_model-class}
\alias{CoModes_model-class}
\title{The [\code{\linkS4class{CoModes_model}}] class}
\description{
  This is a class to contain results from CoModescluster function.
}
\details{
  This class is composed by three attributs
  \describe{
    \item{nbclasses}{the class number of the model.}
    \item{sigma}{this vector indicates the block affectation of the variables.}
    \item{modes}{this matrix indicates the mode number of the multinomial distribution for each block (column) and each class (row) .}
  } 
}
\examples{
  getSlots("CoModes_model")
}
