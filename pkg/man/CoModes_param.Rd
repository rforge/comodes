\name{CoModes_param-class}
\alias{CoModes_param-class}
\title{The [\code{\linkS4class{CoModes_param}}] class}
\description{
  This is a class to contain results from CoModescluster function.
}
\details{
  This class is composed by tow attributs
  \describe{
    \item{proportions}{the vector of the class proportions.}
    \item{alpha}{the list of the parameters of the multinomial. The first level of the list
    corresponds to the class and the second level to the block. Thus, alpha[[ k ]][[ j ]] is a 
    matrix where the first column indicates the mode probability and where the other columns
    indicate the modality crossings for the distribution of the block j for the component k.
    Note that the last row correspond to the probability mass
    uniformly spread among the non-modes modality crossings.}
  } 
}
\examples{
  getSlots("CoModes_param")
}
