\name{SPIAcom}
\alias{SPIAcom}
\title{
Performing topology-based pathway ranking (SPIA) with integrated data
}
\description{
A function performs topology-based pathway ranking (SPIA) with integrated data
}
\usage{
SPIAcom(DEexpression)
}
\arguments{
  \item{DEexpression}{
Six variables including all expressed proteins and phosphoproteins, their differential expression and the fold change of the differential expression.
}
}
\details{
Integrated differential expression is the intersection of the differentially expressed proteins and phosphoproteins. Integrated degree of differential expression is the sum of the fold change absolute value of protein and phosphoptrotein. The integrated information above is the input for SPIA pathway analysis. 2000 bootstrap iterations is used to compute the P PERT value of pathway. Fisher's method is used to combine two kinds of p value from topology of pathway and over representation analysis. The protein/phosphoprotein ID is Entrez gene ID.
}
\value{
The rankings and the significance of pathways with SPIA method.
}
\author{
Jie Ren, Jing Li
}
\seealso{
help(SPIA)
}
\examples{

data(exprs)
data(phosexprs)
case=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10")
ctrl=c("B1","B2","B3")
DEexpression=Dexpress(exprs,phosexprs,case,ctrl)
SPIAresults <- SPIAcom(DEexpression)
head(SPIAresults)
}

