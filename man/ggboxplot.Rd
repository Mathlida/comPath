\name{ggboxplot}
\alias{ggboxplot}
\title{
Boxplot of visualizing different impacts between proteins and phosphoproteins on a disturbed pathway.
}
\description{
A function to plot boxplot which visualzes the different impacts between proteins and phosphoproteins on a disturbed pathway.
}
\usage{
ggboxplot(pathdata,pathname)
}
\arguments{
  \item{pathdata}{
A matrix. Each row represents protein or phosphoprotein included the pathway, and columns represent their fold change of protein and phosphoprotein  respectively. If the protein or phosphoprotedoes not express, the value will be filled with "NA".
}
  \item{pathname}{KEGG ID of the pathtway}
}
\details{
If the impacts between protome and phosphoproteome on a disturbed pathway are different, the function will be called. Users will get The boxplot of the z-scored fold change of protein and phosphoprotein expression included in the pathway. The outliers of proteins or phosphoproteins will be labeled with Entrez gene ID.
}

\author{
Jie Ren, Jing Li
}

\seealso{
help(geom_boxplot)
}
\examples{
ggboxplot(pathdata,pathname)
}

\keyword{ pathdata }
\keyword{ pathname }
