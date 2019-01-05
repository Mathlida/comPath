\name{GeneSetcom}
\alias{GeneSetcom}
\title{
Performing Geneset-based pathway ranking (ORA & GSEA) with integrated data
}
\description{
A function performs Geneset-based pathway ranking (ORA & GSEA) with integrated data
}
\usage{
GeneSetcom(DEexpression)
}
\arguments{
  \item{DEexpression}{
Six variables including all expressed proteins and phosphoproteins, their differential expression and the fold change of the differential expression.
}
}
\details{
Intersection of the differentially expressed proteins and phosphoproteins, the sum of the fold change absolute value of protein and phosphoptrotein are set as the input for Geneset-based pathway analysis. Permutations for deriving p-values in Geneset-based pathway analysis is 2000. The p value adjust method is "Benjamini & Hochberg" (p < 0.05). The protein/phosphoprotein ID is Entrez gene ID.}
\value{
The rankings and the significance of pathways with ORA and GSEA methods.
}
\author{
Jie Ren, Jing Li
}
\seealso{
help(HTSanalyzeR)
}
\examples{
Gsetresults <- GeneSetcom(DEexpression)
ORAresult <- Gsetresults@result$HyperGeo.results$PW_KEGG
GSEAresult <- Gsetresults@result$GSEA.results$PW_KEGG
}

\keyword{ORA}
\keyword{GSEA}
