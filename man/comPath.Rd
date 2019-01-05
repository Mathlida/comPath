\name{comPath-package}
\alias{comPath-package}
\alias{comPath}
\docType{package}
\title{
Integrating proteomic and phosphoproteomic data for pathway analysis.
}
\description{
A package Integrates proteomic and phosphoproteomic data for pathway analysis, evaluates and visualizes the different imapcts between proteome and phosphoproteome on a disturbed pathway.
}
\details{
\tabular{ll}{
Package: \tab SeqMADE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-08-14\cr
LazyLoad: \tab yes\cr
}
The main functions in this package are:
\code{\link{Dexpress}} is a function of identifying the differentially expressed proteins and phosphoproteins and their fold change.
\code{\link{GeneSetcom}} is a function of performing Geneset-based pathway ranking (ORA & GSEA) with integrated data.
\code{\link{SPIAcom}} is a function of performing topology-based pathway ranking (SPIA) with integrated data.
\code{\link{impact}} is a function of evaluating and visualizing the different impacts between proteins and phosphoproteins on a disturbed pathway.
}

\author{
Jie Ren, Jing Li.
Maintainer: JieRen<ntrenjie@163.com>
}
\references{
Ren, J., Wang, B. and Li, J.Integrating proteomic and phosphoproteomic data for pathway analysis in breast cancer. BMC Systems Biology2018,12 (Suppl 8):130
}
\seealso{
help(SPIA), help(HTSanalyzeR)
}
\examples{
data(exprs)
data(phosexprs)
case=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10")
ctrl=c("B1","B2","B3")
DEexpression=Dexpress(exprs,phosexprs,case,ctrl)
summary(DEexpression)
Gsetresults <- GeneSetcom(DEexpression)
ORAresult <- Gsetresults@result$HyperGeo.results$PW_KEGG
GSEAresult <- Gsetresults@result$GSEA.results$PW_KEGG
SPIAresults <- SPIAcom(DEexpression)
}
\keyword{package}

