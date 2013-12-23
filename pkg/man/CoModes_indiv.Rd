\name{CoModes_indiv-class}
\alias{CoModes_indiv-class}
\title{The [\code{\linkS4class{CoModes_indiv}}] class}
\description{
  This is a class to contain results from CoModescluster function.
}
\details{
  This class is composed by three attributs
  \describe{
    \item{proba}{a matrix where the rows correspond to the individuals and where the columns the class. The element of
    the row i and the column k indicates the value of the probability function of the individual i for the component k.}
    \item{tik}{a matrix where the rows correspond to the individuals and where the columns the class. The element of
    the row i and the column k indicates the probability that the individual i arises from the component k.}
    \item{partition}{this vector indicates the class membership of each individual by applying the Maximum A Posterior rule.}
  } 
}
\examples{
  getSlots("CoModes_indiv")
}
