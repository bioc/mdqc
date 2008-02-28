\name{prcomp.robust}
\alias{prcomp.robust}
\title{Principal Components Analysis using Robust Estimators}
\description{
  A function that performs PCA using the robust estimators
  "S-estimator","MCD" and "MVE".
}
\usage{
prcomp.robust(x, robust = c("S-estimator","MCD", "MVE"),
                        nsamp = 10*nrow(x), \dots)
\method{prcomp}{robust}(x, robust = c("S-estimator","MCD", "MVE"),
                        nsamp = 10*nrow(x), \dots)
}

\arguments{
  \item{x}{a matrix. Contains the data to perform PCA on.}
  \item{robust}{The robust estimator to use. One of \code{"S-estimator"},
    \code{"MCD"}, or \code{"MVE"}. The default robust estimator is the S-estimator with 
    25\% breakdown point.}
  \item{nsamp}{The number of subsamples that the robust estimator should
    use. This defaults to 10 times the number of rows in the matrix.}
  \item{\dots}{Further arguments that can be passed to the robust estimator}
}
\details{
   The calculation is done by a singular value decomposition of the
  robust centered and scaled data matrix, not by using
  \code{eigen} on the covariance matrix.  This
  is generally the preferred method for numerical accuracy.  The
  \code{print} method for the these objects prints the results in a nice
  format and the \code{plot} method produces a scree plot. The
  scree plot can be used to determine the number k of principal components preserved in the
  analysis, looking for the \dQuote{elbow} or 
  the first important bend in the line. A biplot can also be generated 
  to represent the values of the first two principal components (PCs) and the
  contribution of each variable to these components in the same plot (see Supplementary
  Material of Cohen Freue et al. (2007)).

}
\value{
  \code{prcomp.robust} returns a list with class \code{"prcomp"}
  containing the following components:
  
  \item{sdev}{the standard deviations of the principal components (i.e.,
    the square roots of the eigenvalues of the covariance
    matrix calculated using the \code{robust} argument, though the calculation is actually done with the singular
    values of the data matrix).}

  \item{rotation}{the matrix of variable loadings (i.e., a matrix whose
    columns contain the eigenvectors).  The function \code{princomp}
    returns this in the element \code{loadings}.}

  \item{x}{the value of the rotated data (the centered and scaled) data multiplied by the \code{rotation} matrix) is returned.}
}

\author{Justin Harrington \email{harringt@stat.ubc.ca} and Gabriela
  V. Cohen Freue \email{gcohen@stat.ubc.ca}.}

\seealso{\code{\link{mdqc}}, \code{prcomp}}
\references{
  Cohen Freue, G. V. and Hollander, Z. and Shen, E. and Zamar, R. H. and Balshaw,
  R. and Scherer, A. and McManus, B. and Keown, P. and McMaster, W. R. and Ng,
  R. T. (2007) \sQuote{MDQC: A New Quality Assessment Method for Microarrays
  Based on Quality Control Reports}. \emph{Bioinformatics} \bold{23}, 3162 -- 3169.
}

\examples{
data(allQC) ## Loads the dataset allQC

prout <- prcomp.robust(allQC)
screeplot(prout, type="line")
biplot(prout)

prout <- prcomp.robust(allQC, robust="MCD")
screeplot(prout, type="line")
biplot(prout)

prout <- prcomp.robust(allQC, robust="MVE")
screeplot(prout, type="line")
biplot(prout)
}
\keyword{multivariate}
\keyword{robust}
