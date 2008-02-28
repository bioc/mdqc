\name{plot.mdqc}
\alias{plot.mdqc}
\title{The plot method for a MDQC object}
\description{
The plot method for a MDQC object, which plots ...
}

\usage{
\method{plot}{mdqc}(x, levels = c(0.9, 0.95, 0.99), xlab="", ylab="",
mfrow=NULL, mfcol=NULL, ...)
}
\arguments{
  \item{x}{An object of the class \sQuote{mdqc}.}
  \item{levels}{A vector or scalar between 0 and 1 for displaying
    critical values for outliers. See details.}
  \item{xlab}{The label for for x-axis. Note that when there are
    multiple plots, the same value of this argument is used for each
    one.}
  \item{ylab}{The label for the y-axis. Note that when there are
    multiple plots, the same value of this argument is used for each
    one.}
  \item{mfrow}{Specify the arrangement of plots on the page, by rows, or
  leave NULL to let the function work it out}
  \item{mfcol}{As for mcol, but arrange plots by column instead}
  \item{\dots}{Other arguments passed to the default plot method.}
}
\details{
  This plot method is for the output from the function
 \code{\link{mdqc}}, and plots the Mahalanobis distances for each
 array. The \code{levels} argument plots horizontal lines at critical
 values (based on the quantiles of a chi-squard distribution), and aids
 in identifying outliers.

For further details, see \cite{Cohen Freue et al. (2007)}

}
\author{Justin Harrington \email{harringt@stat.ubc.ca} and Gabriela
  V. Cohen Freue \email{gcohen@stat.ubc.ca}.}
\references{
 Cohen Freue, G. V. and Hollander, Z. and Shen, E. and Zamar, R. H. and Balshaw,
  R. and Scherer, A. and McManus, B. and Keown, P. and McMaster, W. R. and Ng,
  R. T. (2007) \sQuote{MDQC: A New Quality Assessment Method for Microarrays
  Based on Quality Control Reports}. \emph{Bioinformatics} \bold{23}, 3162 -- 3169.
}

\seealso{\code{\link{mdqc}}}
\examples{
data(allQC)
mdout <- mdqc(allQC, method="cluster", k=3)
plot(mdout)

## Just one critical value
plot(mdout, levels=0.9)
}
\keyword{multivariate}
\keyword{robust}

