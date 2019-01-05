\name{exprs}
\alias{exprs}
\docType{data}
\title{
Protein Expression  Dataset
}
\description{
Protein expression dataset, containing 620 proteins (Entrez gene ID) and 13 samples (10 case and 3 control respectively).
}
\usage{data("exprs")}
\format{
  A data frame with 620 observations on the following 13 variables.
  \describe{
    \item{\code{A1,A2, A3, A4, A5, A6, A7, A8, A9, A10}}{represet case sample names}
    \item{\code{B1, B2, B3}}{represent control sample names}
  }
}
\details{
Protein expression dataset, containing 620 proteins and 13 samples(10 case and 3 control respectively). The protein ID is Entrez gene ID.
}
\examples{
data(exprs)
head(exprs)
}
\keyword{datasets}
