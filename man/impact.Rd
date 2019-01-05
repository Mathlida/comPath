\name{impact}
\alias{impact}
\title{
Evaluate the different impacts betwween proteins and phosphoproteins on a disturbed pathway.
}
\description{
A function evaluates and visualizes the different impacts between proteins and phosphoproteins on a disturbed pathway.
}
\usage{
impact(DEexpression)
}

\arguments{
  \item{DEexpression}{
Six variables including all expressed proteins and phosphoproteins, their differential expression and the fold change of the differential expression.
}
}
\details{
A function evaluates and visualizes the different impacts between proteins and phosphoproteins on a disturbed pathway. The function evaluates the difference by comparing z-scored fold change of proteins and phosphoproteins included in the pathway.If the difference exists in a pathway, z-scored fold change of expressed proteins and phosphoproteins of the pathway will be recorded and a boxplot and a KEGG pathway map will be produced to show the difference. In the pathway mapping figure, the left part of rectangle representing the z-scored fold change of proteins and right part representing phosphoproteins. Red represents up-regulated and yellow represents down-regulated. If the gene does not expressed, the color of the rectangle is white. A file named "plot" will be created in the work directory and all the generated 'png' figures will be included in.
}
\value{
\item{DiffPath}{The number of pathways in which proteins and phosphoproteins contribute differently}
\item{NODiffPath}{The number of pathways in which proteins and phosphoproteins contribute almost equally}
\item{Dlist}{the detail of "DiffPath"}
\item{Slist}{the detail of "NODiffPath"}
}
\author{
Jie Ren, Jing Li
}
\examples{
data(exprs)
data(phosexprs)
case=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10")
ctrl=c("B1","B2","B3")
DEexpression=Dexpress(exprs,phosexprs,case,ctrl)
Impacts <- impact(DEexpression)
summary(Impacts)
}
\keyword{DEexpression}

