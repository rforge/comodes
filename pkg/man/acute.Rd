\encoding{utf8}
\docType{data}
\name{acute}
\alias{acute}
\title{Acute data set}
\format{A data frame with 120 rows and 8 columns.}
\description{
The main idea of this data set is to prepare the algorithm of the expert 
system, which will perform the presumptive diagnosis of two diseases of 
urinary system. It will be the example of diagnosing of the acute 
inflammations of urinary bladder and acute nephritises. For better 
understanding of the problem let us consider definitions of both 
diseases given by medics. Acute inflammation of urinary bladder is
characterised by sudden occurrence of pains in the abdomen region and 
the urination in form of constant urine pushing, micturition pains and 
sometimes lack of urine keeping. Temperature of the body is rising, 
however most often not above 38C. The excreted urine is turbid and 
sometimes bloody. At proper treatment, symptoms decay usually within 
several days. However, there is inclination to returns. At persons with 
acute inflammation of urinary bladder, we should expect that the illness 
will turn into protracted form.

Acute nephritis of renal pelvis origin occurs considerably more often at 
women than at men. It begins with sudden fever, which reaches, and 
sometimes exceeds 40C. The fever is accompanied by shivers and one- or 
both-side lumbar pains, which are sometimes very strong. Symptoms of 
acute inflammation of urinary bladder appear very often. Quite not 
infrequently there are nausea and vomiting and spread pains of whole 
abdomen.

The data was created by a medical expert as a data set to test the 
expert system, which will perform the presumptive diagnosis of two 
diseases of urinary system. The basis for rules detection was Rough Sets 
Theory. Each instance represents an potential patient.

The data was created by a medical expert as a data set to test 
the expert system, which will perform the presumptive diagnosis of two 
diseases of the urinary system.
}
\arguments{
\item{Tem}{Temperature of patient { 35C-42C }}
\item{Nau}{Occurrence of nausea { yes, no }}
\item{Lum}{Lumbar pain { yes, no }} 
\item{Pus}{Urine pushing (continuous need for urination) { yes, no }} 
\item{Mic}{Micturition pains { yes, no }} 
\item{Bur}{Burning of urethra, itch, swelling of urethra outlet { yes, no }} 
\item{Inf}{decision: Inflammation of urinary bladder { yes, no }} 
\item{Nep}{decision: Nephritis of renal pelvis origin { yes, no } } 
}
\examples{
data("acute")
}
\references{

J.Czerniak, H.Zarzycki, Application of rough sets in the presumptive 
diagnosis of urinary system diseases, Artifical Inteligence and Security 
in Computing Systems, ACS'2002 9th International Conference Proceedings, 
Kluwer Academic Publishers,2003, pp. 41-51 
}
\keyword{datasets}
