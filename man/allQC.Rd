\name{allQC}
\alias{allQC}
\docType{data}
\title{QC report for MLL.B}

\description{A subset of arrays from a large acute lymphoblastic leukemia (ALL) study}

\usage{data(allQC)}

\format{
  A data frame with 20 observations on the following 11 variables.
  \describe{
    \item{\code{Scale Factor}}{a numeric vector}
    \item{\code{Percent Present}}{a numeric vector}
    \item{\code{Average Background}}{a numeric vector}
    \item{\code{Minimum Background}}{a numeric vector}
    \item{\code{Maximum Background}}{a numeric vector}
    \item{\code{BioB}}{a numeric vector}
    \item{\code{BioC}}{a numeric vector}
    \item{\code{BioD}}{a numeric vector}
    \item{\code{CreX}}{a numeric vector}
    \item{\code{AFFX-HSAC07/X00351.3'/5'}}{a numeric vector}
    \item{\code{AFFX-HUMGAPDH/M33197.3'/5'}}{a numeric vector}
  }
}

\details{
  Contains the QC report obtained using Bioconductor's simpleaffy
  package for a subset of arrays from a large acute lymphoblastic leukemia
  (ALL) study (\cite{Ross et al., 2004}). The QC report in allQC has
  been generated using the following R commands:

  \preformatted{
    library("affy")

    ## Get the raw data (see help("MLL.B") for further details)
    library("ALLMLL")
    data(MLL.B)

    ## Generate the QC metrics
    library("simpleaffy")
    data.all <- MLL.B[,1:20]
    all.qc <- qc(data.all)

    ## Select relevant information
    allQC <- cbind(sfs(all.qc),percent.present(all.qc)/100,
                   avbg(all.qc),minbg(all.qc),maxbg(all.qc),
                   spikeInProbes(all.qc),ratios(all.qc)[,c(1,3)])
		   
    ## Specify row and column names		   
    colnames(allQC) <- c("Scale Factor","Percent Present", 
	"Average Background", "Minimum Background", "Maximum Background", 
	"BioB", "BioC", "BioD", "CreX", "AFFX-HSAC07/X00351.3'/5'", 
	"AFFX-HUMGAPDH/M33197.3'/5'")
    rownames(allQC) <- 1:20
  }
  Versions 1.16.0 of \sQuote{affy}, 1.2.2 of \sQuote{ALLMLL}, and 2.14.05 of \sQuote{simpleaffy} was
	used.  Part of this dataset has been also studied by \cite{Bolstad et al. (2005)} and \cite{Brettschneider et al. (2007).}
}

\source{
    Ross, M. E. and Zhou, X. and Song, G. and Shurtleff, S. A. and
  Girtman, K. and Williams, W. K. and Liu, H. and Mahfouz, R. and
  Raimondi, S. C. and Lenny, N. and Patel, A. and Downing, J. R. (2003)
  \sQuote{Classification of pediatric acute lymphoblastic leukemia by
    gene expression profiling.} \emph{Blood} \bold{102}, 2951--9.

  All CEL files are freely available in \url{http://www.stjuderesearch.org/data/ALL3/rawFiles.html.}
}
\references{
  Bolstad, B. M. and Collin, F. and Brettschneider, J. and Simpson,
  K. and Cope, L. and Irizarry R. A. and Speed T. P. (2005)
  \sQuote{Quality assessment of Affymetrix GeneChip data.} In Gentleman,
  R. and Carey, C. J. and Huber, W. and Irizarry, R. A. and Dudoit,
  S. \emph{Bioinformatics and Computational Biology Solutions Using R
    and Bioconductor.} New York: Springer.
  
  Brettschneider, J. and Collin, F. and Bolstad, B. M. and Speed,
  T. P. (2007) \sQuote{Quality assessment for short oligonucleotide
    arrays}. Forthcoming in \emph{Technometrics (with Discussion)}.
}

\examples{
data(allQC)
}
\keyword{datasets}
