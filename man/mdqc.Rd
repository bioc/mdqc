\name{mdqc}
\alias{mdqc}
\title{MDQC: Mahalanobis Distance Quality Control}
\description{
  MDQC is a multivariate quality assessment method for
  microarrays based on quality control (QC) reports.
}

\usage{
mdqc(x, method=c("nogroups", "apriori", "global", "cluster", "loading"),
     groups=NULL, k=NULL, pc=NULL,
     robust=c("S-estimator","MCD", "MVE"), nsamp=10*nrow(x))
}

\arguments{
  \item{x}{a numeric matrix or data frame containing the quality
    measures (columns) for each array (rows). The number of rows must
  exceed the number of columns.}
  
  \item{method}{The Mahalanobis Distances (MDs) can be computed on all
    the quality measures in the QC report (this is the default method
    given by \code{method="nogroups"}), on the first k principal components
    resulting from a principal component analysis (PCA) of the QC report
    (\code{"global"}) or on subsets of quality measures in the QC report
    (\code{"apriori"}: groups defined by the user, \code{"cluster"}:
    groups resulting from a cluster analysis, or \code{"loading"}:
    groups resulting from a cluster analysis in the space of the
    loadings of a PCA). While the first two methods compute a single MD
    for each array, the last three compute one MD within each created
    group of quality measures.}
  
  \item{groups}{A list to specify the groups of quality measures when
    the \dQuote{apriori} method is chosen.  E.g. \code{groups =
    list(c(1,2), c(4,6))} puts column 1,2 as one group and 4,6 as a
    second.}
  
  \item{k}{An integer to specify the number of clusters (or groups) to
    be used in the cluster analysis when \dQuote{cluster} or
    \dQuote{loading} methods are chosen.}
  
  \item{pc}{An integer to specify the number of principal components
    analyzed from the PCA when \dQuote{global} or \dQuote{loading}
    methods are chosen.}
  
  \item{robust}{A robust multivariate location/spread estimator (choice
    of S-estimator, MCD or MVE). The default method uses S-estimators
    with a 25\% breakdown point.}

  \item{nsamp}{The number of subsamples that the robust estimator should
    use. This defaults to 10 times the number of rows in the matrix.}
}


\details{
  MDQC flags potentially low quality arrays based on the idea of
  outlier detection, that is, it flags those arrays whose quality
  attributes jointly depart from those of the bulk of the data.
  
  This function computes a distance measure, the Mahalanobis Distance, to
  summarize the quality of each array.  The use of this distance allows us
  to perform a multivariate analysis of the information in QC reports
  taking the correlation structure of the quality measures into
  account. In addition, by using robust estimators to identify the typical
  quality measures of good-quality arrays, the evaluation is not affected
  by the measures of outlying arrays.
  
  MDQC can be based on all the quality measures simultaneously (using
  \code{method="nogroups"}), on subsets of them (using \code{method="apriori"},
  \code{"cluster"}, or \code{"loading"}), or on a transformed space with a lower
  dimension (using \code{method="global"}).
  
  In the \dQuote{apriori} approach the user forms groups of quality
  measures on the basis of an a priori interpretation of them and
  according to the quality aspect they represent.  The \dQuote{cluster}
  and the \dQuote{loading} methods are two data-driven methods to form
  the groups. The former groups the quality measures using clustering
  analysis, and the latter uses the loadings of a principal component
  analysis to identify the quality measures that contain similar
  information and group them. It is important to note that the
  \dQuote{apriori}, the \dQuote{cluster}, and the \dQuote{loading}
  methods create groups of the original quality measures of the report
  and compute one MD within each group. Finally, the \dQuote{global}
  method computes a single MD based on the reduced space of the first k
  principal components from a robust PCA. The number k of PCs can be
  chosen using a scree plot.
  
  More details on each method are given in \cite{Cohen Freue et
    al. (2007)}
}


\value{
  An object of class \sQuote{\dQuote{mdqc}} (with associated plot, print
  and summary methods) with components
  \item{ngroups}{Number of groups in which the MDs have been computed}
  
  \item{groups}{column numbers corresponding to the quality measures
    in each group}
  
  \item{mdqcValues}{Mahalanobis Distance(s) for each array}
  
  \item{x}{dataset containing the numeric quality measures in the
    report}

  \item{method}{method used to group or transform the quality
    measures before computing the MD for each array}
  
  \item{pc}{number of principal components used in the robust PCA.}
  
  \item{k}{number of clusters used in the cluster analysis.}
}

\references{
 Cohen Freue, G. V. and Hollander, Z. and Shen, E. and Zamar, R. H. and Balshaw,
  R. and Scherer, A. and McManus, B. and Keown, P. and McMaster, W. R. and Ng,
  R. T. (2007) \sQuote{MDQC: A New Quality Assessment Method for Microarrays
  Based on Quality Control Reports}. \emph{Bioinformatics} \bold{23}, 3162 -- 3169.

  
  Bolstad, B. M. and Collin, F. and Brettschneider, J. and Simpson,
  K. and Cope, L. and Irizarry R. A. and Speed T. P. (2005)
  \sQuote{Quality assessment of Affymetrix GeneChip data.} In Gentleman
  R. and Carey C. J. and Huber W. and Irizarry R. A. and Dudoit
  S. \emph{Bioinformatics and Computational Biology Solutions Using R
  and Bioconductor}. New York: Springer.
  
  Brettschneider, J. and Collin, F. and Bolstad, B. M. and Speed,
  T. P. (2007) \sQuote{Quality assessment for short oligonucleotide
  arrays}. Forthcoming in \emph{Technometrics (with Discussion)}.
  
  Ross, M. E. and Zhou, X. and Song, G. and Shurtleff, S. A. and
  Girtman, K. and Williams, W. K. and Liu, H. and Mahfouz, R. and
  Raimondi, S. C. and Lenny, N. and Patel, A. and Downing, J. R. (2003)
  \sQuote{Classification of pediatric acute lymphoblastic leukemia by
  gene expression profiling.} \emph{Blood} \bold{102}, 2951--9.
}

\author{Justin Harrington \email{harringt@stat.ubc.ca} and Gabriela
  V. Cohen Freue \email{gcohen@stat.ubc.ca}.}

\note{We thank Christopher Croux for providing us a MATLAB code that
  we translated into R to compute the multivariate S-estimator}


\seealso{\code{\link{prcomp.robust}},\code{\link[cluster]{pam}},
  \code{mahalanobis}, \code{\link{allQC}}}

\examples{

data(allQC)

## Contains the QC report obtained using Bioconductor's simpleaffy package
## for a subset of arrays from a large acute lymphoblastic leukemia (ALL)
## study (Ross et al., 2004).
## This dataset has been also studied by Bolstad et al. (2005) and
## Brettschneider et al. (2007).
## For further information see allQC.


#### No Groups method
# Figure 2 in Cohen Freue et al. (2007):
# Results of MDQC based on all measures of the QC report.

mdout <- mdqc(allQC, method="nogroups")
plot(mdout)
print(mdout)
summary(mdout)

#### A-Priori grouping method
# Figure 3 in Cohen Freue et al. (2007):
# Results of MDQC using the apriori grouping method.

mdout <- mdqc(allQC, method="apriori", groups=list(1:5, 6:9, 10:11))
plot(mdout)



#### Global PCA method
# Figure 4 in Cohen Freue et al.(2007):
# Results of MDQC using the global PCA method.

mdout <- mdqc(allQC, method="global", pc=4)
plot(mdout)



#### Clustering grouping method
# Figure 4 in Supplementary Material of Cohen Freue et al. (2007):
# Results of MDQC using a cluster analysis to form
# 3 groups of quality measures.

mdout <- mdqc(allQC, method="cluster", k=3)
plot(mdout)



#### Loading grouping method
# Figure 4 in Supplementary Material of Cohen Freue et al. (2007):
# Results of MDQC using a cluster analysis on the first
# k=4 loading vectors from a robust PCA to form 3 groups of quality measures.

mdout <- mdqc(allQC, method="loading", k=3, pc=4)
plot(mdout)


### To get the raw MD distances
mdout$mdqcValues

}

\keyword{multivariate}
\keyword{robust}

