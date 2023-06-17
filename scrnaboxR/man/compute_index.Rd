\name{Compute_index}
\alias{Compute_index}
\title{Different indices clustering
}
\description{Compute different indices for different clustering.}
\usage{
compute_index(X,Y)
}
\arguments{
  \item{X}{Lables}
  \item{Y}{Labels}
}
\details{
This function compute indices for the given clustersing; the output is the contigency table, pearson X-squared, Cramer's V, ACC (Computes one subtract the error rate. It uses the Hungarian method to match clusters),Adjusted Mutual Information (AMI), computes the normalized mutual information (NMI), Rand Index (RI), Adjusted Rand Index (ARI)
}
\author{Saeid Amiri}
\examples{
# X<-sample(c('a1', 'a2'), 1000, replace=TRUE)
# Y<-sample(c('b1', 'b2', 'b3'), 1000, replace=TRUE)
# aa<-compute_index(X,Y)
}
\keyword{Annotation}
