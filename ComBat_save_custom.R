# Modified versions of ComBat function from (SVA) R package
# Allows for extraction and application of custom correction factors
# 
# Modified from Surrogate Variable Analysis (SVA) R package
# Modified by Finsam Samson
#
#
#
#
# Examples of the new functions being used are in lines 463-476 (ComBat_save)
# and lines 778-797 (ComBat_custom)
#
#
#
# ComBat_save runs ComBat, and outputs a corrected counts table,
# gamma correction factor, and delta correction factor in one object
#
# ComBat_custom uses the custom gamma and delta correction factors generated
# from ComBat_save. It outputs a corrected counts table.
#
# The two input datasets for ComBat_save and ComBat_custom must have the same
# number of batches
#
# The other functions are helper functions
#
#


library(BatchQC)
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggplot2)

################################################################################################
################################################################################################
################################################################################################

# Helper functions for ComBat_save and ComBat_custom to work properly

library(BiocParallel) 

sva.class2Model <- function(classes) {
  return(model.matrix(~factor(classes)))
}


modefunc <- function(x) {
  return(as.numeric(names(sort(-table(x)))[1]))
}

mono <- function(lfdr){
  .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
}

edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  
  n <- length(p)
  transf <- match.arg(transf)
  
  if(transf=="probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1-eps)
    x <- qnorm(p)
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0*dnorm(x)/y
  }
  
  if(transf=="logit") {
    x <- log((p+eps)/(1-p+eps))
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1+exp(x))^2
    lfdr <- pi0 * dx/y
  }
  
  if(trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if(monotone) {	
    lfdr <- lfdr[order(p)]
    lfdr <- mono(lfdr)
    lfdr <- lfdr[rank(p)]
  }
  
  return(lfdr)
}



# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
  tmp <- strsplit(colnames(dat),'\\.')
  tr <- NULL
  for (i in 1:length(tmp)) {
    tr <- c(tr, tmp[[i]][1] != 'X')
  }
  tr
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2*s2 + m^2) / s2
}

bprior <- function(gamma.hat){
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m*s2 + m^3) / s2
}

postmean <- function(g.hat,g.bar,n,d.star,t2){
  (t2*n*g.hat + d.star*g.bar) / (t2*n + d.star)
}

postvar <- function(sum2,n,a,b){
  (.5*sum2 + b) / (n/2 + a - 1)
}

# Inverse gamma distribution density function. (Note: does not do any bounds checking on arguments)
dinvgamma <- function (x, shape, rate = 1/scale, scale = 1) {
  # PDF taken from https://en.wikipedia.org/wiki/Inverse-gamma_distribution
  # Note: alpha = shape, beta = rate
  stopifnot(shape > 0)
  stopifnot(rate > 0)
  ifelse(x <= 0, 0, ((rate ^ shape) / gamma(shape)) * x ^ (-shape - 1) * exp(-rate/x))
}

# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1,ncol(sdat))))^2, na.rm=TRUE)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  ## cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

## likelihood function used below
L <- function(x,g.hat,d.hat){
  prod(dnorm(x, g.hat, sqrt(d.hat)))
}

## Monte Carlo integration functions
int.eprior <- function(sdat, g.hat, d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2 %*% j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star, sum(g*LH)/sum(LH))
    d.star <- c(d.star, sum(d*LH)/sum(LH))
    ## if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust	
} 

## fits the L/S model in the presence of missing data values

Beta.NA <- function(y,X){
  des <- X[!is.na(y),]
  y1 <- y[!is.na(y)]
  B <- solve(crossprod(des), crossprod(des, y1))
  B
}

################################################################################################
################################################################################################
################################################################################################

ComBat_save <- function(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                   mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) {
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!  
  
  ## coerce dat into a matrix
  dat <- as.matrix(dat)
  
  ## find genes with zero variance in any of the batches
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level){
    if(sum(batch==batch_level)>1){
      return(which(apply(dat[, batch==batch_level], 1, function(x){var(x)==0})))
    }else{
      return(which(rep(1,3)==2))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", length(zero.rows)))
    # keep a copy of the original data matrix and remove zero var rows
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
  
  ## make batch a factor and make a set of indicators for batch
  if(any(table(batch)==1)){mean.only=TRUE}
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch =",ref.batch, "as a reference batch (this batch won't change)")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  message('Standardizing Data across genes')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  
  ##Find Priors
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, aprior) # FIXME 
  b.prior <- apply(delta.hat, 1, bprior) # FIXME
  
  ## Plot empirical and parametric priors
  
  if (prior.plots && par.prior) {
    old_pars <- par(no.readonly = TRUE)
    on.exit(par(old_pars))
    par(mfrow=c(2,2))
    
    ## Top left
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    
    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)
    
    ## Bottom Left
    tmp <- density(delta.hat[1,])
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)
    
    ## Bottom Right
    invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }
  
  ## Find EB batch adjustments
  
  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                       delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                       b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star=temp[1,], delta.star=temp[2,])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    ## DEBUG LINES
    samples_umrr <- bayesdata[,i]
    sub_umrr <- t(batch.design[i,]%*%gamma.star)
    ##
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  ## put genes with 0 variance in any batch back in data
  if (length(zero.rows) > 0) {
    dat.orig[keep.rows, ] <- bayesdata
    bayesdata <- dat.orig
  }
  
  ## EXTRACT
  gamma.star.out <- gamma.star
  delta.star.out <- delta.star
  colnames(gamma.star.out) <- rownames(s.data)
  colnames(delta.star.out) <- rownames(s.data)
  combat_save_out <- list("data"=bayesdata, "gamma"=gamma.star.out, "delta"=delta.star.out)
  
  # return(bayesdata)
  return(combat_save_out)
}

################################################################################################
################################################################################################
################################################################################################

# Example usage of ComBat_save

# Import unnormalized UMRR counts as data frames,
# converts to metrices, and merges them
# raw_242_umrr <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-242-245_FLT-GC/URR/GLDS-242_UMRR_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
# raw_245_umrr <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-242-245_FLT-GC/URR/GLDS-245_UMRR_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
# 
# raw_umrr <- cbind(raw_242_umrr, raw_245_umrr)
# 
# bc_combat <- c(242,242,242,242,242,242,242,242,242,245,245,245,245,245,245)
# combat_umrr_all <- ComBat_save(dat=raw_umrr, batch=bc_combat, par.prior=TRUE, prior.plots=FALSE)
# combat_umrr <- combat_umrr_all$data
# umrr_gamma <- combat_umrr_all$gamma
# umrr_delta <- combat_umrr_all$delta

################################################################################################
################################################################################################
################################################################################################

ComBat_custom <- function(dat, batch, custom.gamma, custom.delta, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                        mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) {
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!  
  
  ## coerce dat into a matrix
  dat <- as.matrix(dat)
  
  ## find genes with zero variance in any of the batches
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level){
    if(sum(batch==batch_level)>1){
      return(which(apply(dat[, batch==batch_level], 1, function(x){var(x)==0})))
    }else{
      return(which(rep(1,3)==2))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", length(zero.rows)))
    # keep a copy of the original data matrix and remove zero var rows
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
  
  ## make batch a factor and make a set of indicators for batch
  if(any(table(batch)==1)){mean.only=TRUE}
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch =",ref.batch, "as a reference batch (this batch won't change)")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  message('Standardizing Data across genes')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  
  ##Find Priors
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, aprior) # FIXME 
  b.prior <- apply(delta.hat, 1, bprior) # FIXME
  
  ## Plot empirical and parametric priors
  if (prior.plots && par.prior) {
    old_pars <- par(no.readonly = TRUE)
    on.exit(par(old_pars))
    par(mfrow=c(2,2))
    
    ## Top left
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    
    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)
    
    ## Bottom Left
    tmp <- density(delta.hat[1,])
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)
    
    ## Bottom Right
    invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }
  
  ## Find EB batch adjustments
  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                       delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                       b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star=temp[1,], delta.star=temp[2,])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  
  
  ## CREATE BLANKS
  ## gamma is additive, delta is multiplicative
  gamma.star.blank <- gamma.star
  colnames(gamma.star.blank) <- rownames(s.data)
  gamma.star.blank[] <- 0
  delta.star.blank <- delta.star
  colnames(delta.star.blank) <- rownames(s.data)
  delta.star.blank[] <- 1
  
  ## Replace matching columns
  # Select indices, omit nans
  gamma_indices1_with_nan <- match(colnames(custom.gamma), colnames(gamma.star.blank))
  gamma_indices1 <- na.omit(gamma_indices1_with_nan)
  gamma_indices2_with_nan <- match(colnames(gamma.star.blank), colnames(custom.gamma)) # GETS NEEDED INDICES TO INDEX FROM CUSTOM
  gamma_indices2 <- na.omit(gamma_indices2_with_nan)
  
  delta_indices1_with_nan <- match(colnames(custom.delta), colnames(delta.star.blank))
  delta_indices1 <- na.omit(delta_indices1_with_nan)
  delta_indices2_with_nan <- match(colnames(delta.star.blank), colnames(custom.delta)) # GETS NEEDED INDICES TO INDEX FROM CUSTOM
  delta_indices2 <- na.omit(delta_indices2_with_nan)
  
  # Replace
  gamma.star.blank[,gamma_indices1] <- custom.gamma[,gamma_indices2]
  delta.star.blank[,delta_indices1] <- custom.delta[,delta_indices2]
  
  j <- 1
  for (i in batches){
    ## DEBUG LINES
    #samples_run <- bayesdata[,i]
    #sub_run <- t(batch.design[i,]%*%gamma.star)
    ##
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star.blank))/(sqrt(delta.star.blank[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  ## put genes with 0 variance in any batch back in data
  if (length(zero.rows) > 0) {
    dat.orig[keep.rows, ] <- bayesdata
    bayesdata <- dat.orig
  }
  
  # combat_save_out <- list("data"=bayesdata, "gamma"=gamma.star, "delta"=delta.star)
  
  return(bayesdata)
  # return(combat_save_out)
}

################################################################################################
################################################################################################
################################################################################################

# Example usage of ComBat_custom

# Import unnormalized raw counts as data frames
# Extract flight and ground control samples, combine
# raw_242 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-242/GLDS-242_input_files/GLDS-242_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# raw_242_flt <- as.matrix(raw_242[,grep("FLT", names(raw_242))])
# raw_242_gc <- as.matrix(raw_242[,grepl("GC", names(raw_242))])
# raw_242_flt_gc <- cbind(raw_242_flt, raw_242_gc)
# 
# raw_245 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-245/GLDS-245_input_files/GLDS-245_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# raw_245_flt <- as.matrix(raw_245[,grepl("FLT", names(raw_245))])
# raw_245_gc <- as.matrix(raw_245[,grepl("GC", names(raw_245))])
# raw_245_flt_gc <- cbind(raw_245_flt, raw_245_gc)
# 
# raw_counts <- cbind(raw_242_flt_gc, raw_245_flt_gc)
# 
# bc_combat <- c(rep.int(242, 9),rep.int(245, 39))
# cov_combat <- c(rep.int(1, 5),rep.int(2, 4),rep.int(1, 20),rep.int(2, 19))
# combat_man_counts <- ComBat_custom(dat=raw_counts, batch=bc_combat, custom.gamma=umrr_gamma, custom.delta=umrr_delta, mod=cov_combat, par.prior=TRUE, prior.plots=FALSE)
# combat_man_counts is just the corrected counts table, just like what the output of ComBat is

################################################################################################
################################################################################################
################################################################################################

