% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DirichletMultinomialClassifier.R
\name{DirichletMultinomialClassifier}
\alias{DirichletMultinomialClassifier}
\title{Dirichlet multinomial classifier}
\usage{
DirichletMultinomialClassifier(X, catalogue.profile, catalogue.celltype,
  catalogue.prior = NULL, X.cluster = NULL, genes.use = NULL,
  cells.use = NULL, do.training = T, BPPARAM = SerialParam())
}
\arguments{
\item{BPPARAM}{}
}
\value{

}
\description{
\code{DirichletMultinomialClassifier} computes the predictive probability for the given cells following
Multinomial distribution model with the Dirichlet conjugate prior.
}
\details{

}
\examples{

}
\author{
Feng Zeng <zengfeng@xmu.edu.cn>
}
