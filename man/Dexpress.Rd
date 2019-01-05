\name{Dexpress}
\alias{Dexpress}
\title{
Identifying the differentially expressed proteins and phosphoproteins and the differentially expressed degree (fold change) of all the proteins and phosphoproteins
}
\description{
Identifying the differentially expressed proteins and phosphoproteins and the differentially expressed degree (fold change) of all the proteins and phosphoproteins
}
\usage{
Dexpress(exprs,phosexprs,case,ctrl)
}
\arguments{
  \item{exprs}{The protein expression values which have been .}
  \item{phosexprs}{The phosphoproteinexpression values}
  \item{case}{A vector of names of case samples}
  \item{ctrl}{A vector of names of ctrl samples}
}
\details{
If this function is called, two queries will appear in the console: "Have you normalized exprs dataset? Yes(y) or NO(n)?" & "Have you normalized phosexprs dataset? Yes(y) or NO(n)?". If you type "y", Student’s t-test will be used to get the differential expression of proteome or phosphoproteome. If you type "n",  ‘Mann-Whitney’ test will be used to get the differential expression. A Bonferroni correction (p < 0.05) will be applied following the two kinds of tests. The differentially expressed degree (fold change) of all the proteins and phosphoproteins is between case and control samples. The input of the function are  proteome and phosphoproteome expression data matrix and two vector of sample names (case & control).
}
\value{
\item{allpro}{all of the expressed proteins}
\item{allphos}{all of the expressed phosphoproteins}
\item{DEpro}{Differentilly expressed proteins}
\item{DEphos}{Differentilly expressed phosphoproteins}
\item{fcpro}{the fold chage of DEpros}
\item{fcphos}{the fold chage of DEphos}
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
summary(DEexpression)
}
\keyword{DEpro}
\keyword{DEphos}
\keyword{fcpro}
\keyword{fcall_phos}
\keyword{allpro}
\keyword{allphos}
