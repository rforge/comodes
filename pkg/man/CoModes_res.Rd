\name{CoModes_res-class}
\alias{CoModes_res-class}
\title{The [\code{\linkS4class{CoModes_res}}] class}
\description{
  This is a class to contain results from CoModescluster function.
}
\details{
  This class is composed by five attributs
  \describe{
  \item{data}{an instance of the [\code{\linkS4class{CoModes_data}}]  class containing the elements relative to the data.}
  \item{model}{an instance of the [\code{\linkS4class{CoModes_model}}]  class containing the elements relative to the model.}
  \item{param}{an instance of the [\code{\linkS4class{CoModes_param}}]  class containing the elements relative to the parameters.}
  \item{indiv}{an instance of the [\code{\linkS4class{CoModes_indiv}}]  class containing the elements relative to the individuals.}
  \item{criteria}{an instance of the [\code{\linkS4class{CoModes_criteria}}]  class containing the elements relative to the criteria.}
    } 
}
\examples{
getSlots("CoModes_res")
}
