##  code by Justin Harrington

## mdqc - Methods:
##  "nogroups" - use all columns
##  "apriori" - use groups of columns (elements of a list)
##  "global" - global pca (need to specify how many principal components)
##  "cluster" - use (pam) to cluster columns (need to specify how many clusters)
##  "loading" - loading PCA method by clustering in the space of the loadings (need to specify how many clusters & PCs)

################
## Main Function
################
mdqc <- function(x, method=c("nogroups", "apriori", "global", "cluster", "loading"),
                 groups=NULL, k=NULL, pc=NULL, robust=c("S-estimator","MCD", "MVE"),
                 nsamp=10*nrow(x)){
  ## x - dataset. Matrix or dataframe. Must be numeric
  ## method (see above)
  ## groups - for apriori method, specify groups as a list. E.g. groups = list(c(1,2), c(4,6)) puts column 1,2 as one group and 4,6 as a second
  ## k - number of cluster (cluster method & loading method)
  ## pc - number of principal components to use (global PCA method $ loading method)
  ## robust - which robust measure of location/spread (choice of S-estimator or MCD)
  
  ## Perform argument match for method & rmethod - if doesn't match, will fail
  method <- match.arg(method)
  robust <- match.arg(robust)
  rmethod <- switch(robust, "S-estimator"="mySm",
                    "MCD"="cov.mcd",
                    "MVE"="cov.mve")
  
  ## data validation on x
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("The matrix x does not appear to be numeric.")
  if (dim(x)[1] < dim(x)[2] & method == "nogroups") stop("The matrix x has less rows than columns.")

  ## Argument validation (on k, pc, groups etc)
  if (method=="apriori"){
    if (is.null(groups)) stop("You have selected the a-priori method without specifying groupings.")
    if (any(sapply(groups, length) > nrow(x))) stop("The matrix x has less rows than columns in at least one group.")
  }
  if (method %in% c("global","loading")){
    if (pc > nrow(x)) stop("The matrix x has less rows than principal components.")
    if (is.null(pc)) stop("You need to specify the number of principal components for the ", method, " method.")
    if (pc < 1) stop("The number of principal components must be an integer greater than 0.")
  }
  if (method %in% c("cluster","loading")){
    if (is.null(k)) stop("You need to specifiy the number of clusters for the ", method, " method.")
    if (k < 1) stop("The number of clusters must be an integer greater than 0.")
  }
  
  ## By default
  y <- x 
  ## otherwise create transformed space (in case of global or loading)
  if (method %in% c("global", "loading")){
    ## Perform Robust PCA
    robPCA <- prcomp.robust(x, robust, nsamp=nsamp)
    if (method == "global") ## use this dataset for glocal PCA method
      y <- robPCA$x
  }

  ## create groupings when method is cluster/loading
  if (method == "cluster"){
    robCov <- do.call(rmethod, list(x=y, nsamp=nsamp))$cov
    robCor <- diag(1/sqrt(diag(robCov)))%*%robCov%*%diag(1/sqrt(diag(robCov)))
    robDist <- as.dist((1-robCor)/2)
    clustout <- pam(x=robDist, k=k, cluster.only=TRUE)
    if (any(table(clustout) > nrow(x))) stop("The matrix x has less rows than columns in at least one of the groups from clustering.")
  }

  if (method == "loading"){
    clustout <- pam(x=robPCA$rotation[, (1:pc)], metric="manhattan", k=k, cluster.only=TRUE) 
  }

  ## apply method to form groups
  groups <- switch(method,
                   "nogroups"=list(1:ncol(y)),
                   "global"=list(1:pc),
                   "apriori"=groups,
                   "cluster"=,"loading"=split(1:ncol(y), clustout))
  
  ngroups <- length(groups)

  ## Calculate mahalanobis distances
  mdqcValues <- list()
  for (i in 1:ngroups){
    if (method == "global"){
      robMdqc <- do.call(rmethod, list(x=y, nsamp=nsamp))
      robMdqc$center <- robMdqc$center[1:pc]
      robMdqc$cov <- robMdqc$cov[1:pc, 1:pc]
    }
    else
      robMdqc <- do.call(rmethod, list(x=y[, groups[[i]] ], nsamp=nsamp))
      
    if (length(groups[[i]]) > 1)
      mdqcValues[[i]] <- sqrt(mahalanobis(y[, groups[[i]] ], robMdqc$center, robMdqc$cov))
    else ## If number in group is one
      mdqcValues[[i]] <- abs(scale(y[, groups[[i]] ], center=robMdqc$center, scale=sqrt(robMdqc$cov)))
  }

  ## Return data
  mdqcObject <- list(ngroups=ngroups, groups=groups, mdqcValues=mdqcValues, x=x, method=method, robust=robust, pc=pc, k=k)
  class(mdqcObject) <- "mdqc"
  return(mdqcObject)
}


#####################
## robust prcomp
#####################

prcomp.robust <- function(x, robust=c("S-estimator","MCD", "MVE"), nsamp=10*nrow(x), ...){
  robust <- match.arg(robust)
  rmethod <- switch(robust,
                    "S-estimator"="mySm",
                    "MCD"="cov.mcd",
                    "MVE"="cov.mve") 
  
  ## Standardize data using robust esimates
  robEstimates <- do.call(rmethod, list(x=x, nsamp=nsamp))
  scaledData <- scale(x, center=robEstimates$center, scale=sqrt(diag(robEstimates$cov)))
  
  ## Use robust estimates for calculating the PCs
  robEstimates <- do.call(rmethod, list(x=scaledData, nsamp=nsamp))
  eigenData <- svd(robEstimates$cov)
  y <- scaledData %*% eigenData$v
  outputs <- list(sdev=sqrt(eigenData$d), rotation=eigenData$v, x=y)
  dimnames(outputs$rotation) <- list(colnames(x), paste("PC",1:ncol(x), sep=""))
  class(outputs) <- "prcomp"
  return(outputs)
}


#####################
########## S3 Methods
#####################
plot.mdqc <- function(x, levels=c(0.9, 0.95, 0.99), xlab="", ylab="", mfrow=NULL, mfcol=NULL, ...){
  if (x$ngroups > 1){
    ## Set up mfrow (for multiple plots)
    if (is.null(mfrow) & is.null(mfcol)){
      mfdim <- ceiling(sqrt(x$ngroups)) ## Lazy default
      opar <- par(mfrow=c(ceiling(x$ngroups/mfdim), mfdim))
    }
    else
      opar <- par(mfrow=mfrow, mfcol=mfcol)
    on.exit(par(opar))
  }

  index <- 1:nrow(x$x)

  for (i in 1:x$ngroups){
    r <- length(x$groups[[i]])
    plot(index, x$mdqcValues[[i]], xlab=xlab, ylab=ylab, main=ifelse(x$ngroups > 1, paste("Group",i), ""), ...)
    abline(h=sqrt(qchisq(levels, r)), lwd=1.3, lty=1:3)
    flag <- which(x$mdqcValues[[i]] > sqrt(qchisq(min(levels), r)))
    points(index[flag], x$mdqcValues[[i]][flag], pch=19, col=1)
    text(index[flag], x$mdqcValues[[i]][flag]-par("cxy")[2], flag)
  }
}

summary.mdqc <- function(object, ...) {
  cat("\nSummary information for MDQC\n")
  cat("Method used:", object$method, "\tNumber of groups:", object$ngroups, "\nRobust estimator:", object$robust)
  if (object$method %in% c("loading", "global"))
    cat("\tNumber of Principal Components:", object$pc,"\n")
  cat("\nNumber of Outliers:\n")
  for (i in 1:object$ngroups) {
    if (object$ngroups > 1) {
      cat("Group", i, " - Columns\n")
      print(object$group[[i]])
    }
    r <- length(object$groups[[i]])
    cat("90%\t95%\t99%\n")
    for (j in 1:3){
      crit <- c(0.9, 0.95, 0.99)[j]
      flag <- which(object$mdqcValues[[i]] > sqrt(qchisq(crit, r)))
      cat(length(flag), "\t")
    }
    cat("\n")
  }
}

print.mdqc <- function(x, ... ){
  index <- 1:nrow(x$x)
  cat("Method used:", x$method, "\tNumber of groups:", x$ngroups, "\nRobust estimator:", x$robust)
  if (x$method %in% c("loading", "global"))
    cat("\tNumber of Principal Components:", x$pc,"\n")

  for (i in 1:x$ngroups) {
    if (!(x$method %in% c("global", "nogroups"))) {
      cat("Group", i," - Columns\n")
      print(x$group[[i]])
    }
    r <- length(x$groups[[i]])
    for (j in 1:3){
      crit <- c(0.9, 0.95, 0.99)[j]
      cat("MDs exceeding the square root of the ", crit*100, "% percentile of the Chi-Square distribution\n")
      flag <- which(x$mdqcValues[[i]] > sqrt(qchisq(crit, r)))
      print(index[flag])
    }
    cat("\n")
  }
}

