\name{phosexprs}
\alias{phosexprs}
\docType{data}
\title{
Phosphoprotein Expression  Dataset
}
\description{
Phosphoprotein expression dataset, containing 620 phosphoproteins (Entrez gene ID) and 13 samples(10 case and 3 control respectively).
}
\usage{data("phosexprs")}
\format{
  A data frame with 620 observations on the following 13 variables.
  \describe{
    \item{\code{A1,A2, A3, A4, A5, A6, A7, A8, A9, A10}}{represent case sample names}
    \item{\code{B1, B2, B3}}{represent control sample names}
  }
}
\details{
Phosphoprotein expression dataset, containing 620 phosphoproteins and 13 samples (10 case and 3 control respectively).The phosphoprotein ID is Entrez gene ID.
}
\examples{
data(phosexprs)
head(phosexprs)
}
\keyword{datasets}
