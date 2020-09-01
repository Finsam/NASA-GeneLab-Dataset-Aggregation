# Modified versions of ComBatseq_ function from (SVA) R package
# Allows for extraction and application of custom correction factors
# 
# Modified from Surrogate Variable Analysis (SVA) R package
# Modified by Finsam Samson
#
#
#
#
# Examples of the new functions being used are in lines 331-353 (ComBat_seq_save)
# and lines 575-647 (ComBat_seq_custom)
#
#
#
# ComBat_seq_save runs ComBat_seq, and outputs a corrected counts table,
# gamma correction factor, and phi correction factor in one object
#
# ComBat_seq_custom uses the custom gamma and phi correction factors generated
# from ComBat_seq_save. It outputs a corrected counts table.
#
# The two input datasets for ComBat_seq_save and ComBat_seq_custom must have the same
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
library(edgeR)
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(BiocParallel)

################################################################################################
################################################################################################
################################################################################################

# Helper functions for ComBat_seq_save and ComBat_seq_custom to work properly

####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})  
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }
    
    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)	
  return(res)
} 


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(counts_sub[a, b] <= 1){
        new_counts_sub[a,b] <- counts_sub[a, b]
      }else{
        tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}

################################################################################################
################################################################################################
################################################################################################

ComBat_seq_save <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE, shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){
  ########  Preparation  ########  
  ## Does not support 1 sample per batch yet
  batch <- as.factor(batch)
  if(any(table(batch)<=1)){
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  
  ## Remove genes with only 0 counts in any batch
  keep_lst <- lapply(levels(batch), function(b){
    which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  
  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=counts)
  
  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')
  
  ## Make design matrix 
  # batch
  batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(counts)))
  }
  # drop intercept in covariate model
  if(!is.null(covar_mod)){
    if(is.data.frame(covar_mod)){
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
  }
  # bind with biological condition of interest
  mod <- cbind(mod, covar_mod)
  # combine
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
  }
  
  ## Check for missing values in count matrix
  NAs = any(is.na(counts))
  if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}
  
  
  ########  Estimate gene-wise dispersions within each batch  ########
  cat("Estimating dispersions\n")
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
      # not enough residual degree of freedom
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
    }else{
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
    }
  })
  
  ## Estimate gene-wise dispersion within each batch 
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
      # not enough residual degrees of freedom - use the common dispersion
      return(rep(disp_common[j], nrow(counts)))
    }else{
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                    dispersion=disp_common[j], prior.df=0))
    }
  })
  names(genewise_disp_lst) <- paste0('batch', levels(batch))
  
  ## replace genewise_disp_list
  
  
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
  }
  
  ## doubt about the new offset
  
  ########  Estimate parameters from NB GLM  ########
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
  glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 
  
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  ## DEBUG LINES
  test1 <- glm_f2$coefficients
  
  
  ##
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########  
  if(shrink){
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii){
      if(ii==1){
        mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                           gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
      }else{
        invisible(capture.output(mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                                                    gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))
    
    gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
    phi_star_mat <- do.call(cbind, phi_star_mat)
    
    if(!shrink.disp){
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }else{
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(jj in 1:n_batch){
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], n_batches[jj])) # calculate with new gamma
  }
  phi_star <- rowMeans(phi_star_mat)
  
  ########  Adjust the data  ########  
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts)) # dimensions of cut output (must adjust again)
  for(kk in 1:n_batch){ # keep n_batch the same in new
    counts_sub <- counts[, batches_ind[[kk]]] # counts will mismatch
    old_mu <- mu_hat[, batches_ind[[kk]]] # mu_hat will mismatch (mu may need to be kept)
    old_phi <- phi_hat[, kk] # phi_hat will mismatch
    new_mu <- mu_star[, batches_ind[[kk]]] # mu_star will mismatch (mu may need to be kept)
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub, 
                                                          old_mu=old_mu, old_phi=old_phi, 
                                                          new_mu=new_mu, new_phi=new_phi)
  }
  
  #dimnames(adjust_counts) <- dimnames(counts)
  #return(adjust_counts)
  
  ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
  adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[keep, ] <- adjust_counts
  adjust_counts_whole[rm, ] <- countsOri[rm, ]
  
  ## EXTRACT
  combat_seq_save_out <- list("data"=adjust_counts_whole, "gamma"=gamma_star_mat, "phi"=phi_hat)
  
  # return(adjust_counts_whole)
  return(combat_seq_save_out)
}

################################################################################################
################################################################################################
################################################################################################

# Example usage of ComBat_seq_save

# # Import unnormalized UMRR counts as data frames,
# # converts to metrices, and merges them
# raw_242_umrr <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-242/GLDS-242_URR/GLDS-242_UMRR_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
# raw_245_umrr <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-245/GLDS-245_URR/GLDS-245_UMRR_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
# 
# raw_umrr <- cbind(raw_242_umrr, raw_245_umrr)
# 
# # Create batch covariate table (TODO: make this better - for other datasets
# bc_batchqc <- matrix(c(242,242,242,242,242,242,242,242,242,245,245,245,245,245,245))
# 
# 
# # batchQC(raw_umrr, batch=bc_batchqc, report_file = "batchqc_umrr_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# 
# bc_combat <- c(242,242,242,242,242,242,242,242,242,245,245,245,245,245,245)
# 
# # test_c_seq <- ComBat_seq(counts=raw_umrr, batch=bc_combat)
# combat_seq_umrr_all <- ComBat_seq_save(counts=raw_umrr, batch=bc_combat)
# 
# 
# combat_seq_umrr <- combat_seq_umrr_all$data
# umrr_gamma <- combat_seq_umrr_all$gamma
# umrr_phi <- combat_seq_umrr_all$phi
## batchQC(combat_seq_umrr, batch=bc_batchqc, report_file = "batchqc_umrr_combat_seq_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")

################################################################################################
################################################################################################
################################################################################################

ComBat_seq_custom <- function(counts, batch, custom.gamma, custom.phi, group=NULL, covar_mod=NULL, full_mod=TRUE, shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){  
  ########  Preparation  ########  
  ## Does not support 1 sample per batch yet
  batch <- as.factor(batch)
  if(any(table(batch)<=1)){
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  
  ## Remove genes with only 0 counts in any batch
  keep_lst <- lapply(levels(batch), function(b){
    which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  
  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=counts)
  
  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')
  
  ## Make design matrix 
  # batch
  batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(counts)))
  }
  # drop intercept in covariate model
  if(!is.null(covar_mod)){
    if(is.data.frame(covar_mod)){
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
  }
  # bind with biological condition of interest
  mod <- cbind(mod, covar_mod)
  # combine
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
  }
  
  ## Check for missing values in count matrix
  NAs = any(is.na(counts))
  if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}
  
  
  ########  Estimate gene-wise dispersions within each batch  ########
  cat("Estimating dispersions\n")
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
      # not enough residual degree of freedom
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
    }else{
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
    }
  })
  
  
  ## Estimate gene-wise dispersion within each batch 
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
      # not enough residual degrees of freedom - use the common dispersion
      return(rep(disp_common[j], nrow(counts)))
    }else{
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                    dispersion=disp_common[j], prior.df=0))
    }
  })
  names(genewise_disp_lst) <- paste0('batch', levels(batch))
  
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
  }
  
  
  ########  Estimate parameters from NB GLM  ########
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
  glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 
  
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########  
  if(shrink){
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii){
      if(ii==1){
        mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                           gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
      }else{
        invisible(capture.output(mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                                                    gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))
    
    gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
    phi_star_mat <- do.call(cbind, phi_star_mat)
    
    if(!shrink.disp){
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }else{
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  
  # REPLACE
  # Create blanks (size of current)
  gamma.star.blank <- gamma_star_mat
  rownames(gamma.star.blank) <- rownames(gamma_star_mat)
  colnames(gamma.star.blank) <- colnames(gamma_star_mat)
  gamma.star.blank[] <- 0
  phi.hat.blank <- phi_hat
  rownames(phi.hat.blank) <- rownames(phi_hat)
  colnames(phi.hat.blank) <- colnames(phi_hat)
  phi.hat.blank[] <- 0
  
  # Replace matching columns
  # Select indices, omit nans
  gamma_indices1_with_nan <- match(rownames(custom.gamma), rownames(gamma.star.blank))
  gamma_indices1 <- na.omit(gamma_indices1_with_nan)
  gamma_indices2_with_nan <- match(rownames(gamma.star.blank), rownames(custom.gamma)) # GETS NEEDED INDICES TO INDEX FROM CUSTOM
  gamma_indices2 <- na.omit(gamma_indices2_with_nan)
  
  phi_indices1_with_nan <- match(colnames(custom.phi), colnames(phi.hat.blank))
  phi_indices1 <- na.omit(phi_indices1_with_nan)
  phi_indices2_with_nan <- match(colnames(phi.hat.blank), colnames(custom.phi)) # GETS NEEDED INDICES TO INDEX FROM CUSTOM
  phi_indices2 <- na.omit(phi_indices2_with_nan)
  
  # Replace
  gamma.star.blank[gamma_indices1,] <- custom.gamma[gamma_indices2,]
  phi.hat.blank[phi_indices1,] <- custom.phi[phi_indices2,]
  
  phi_star_mat <- phi.hat.blank
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(jj in 1:n_batch){
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma.star.blank[, jj], n_batches[jj])) # used to be gamma_star_mat
  }
  phi_star <- rowMeans(phi_star_mat)
  
  ########  Adjust the data  ########  
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(kk in 1:n_batch){
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub, 
                                                          old_mu=old_mu, old_phi=old_phi, 
                                                          new_mu=new_mu, new_phi=new_phi)
  }
  
  #dimnames(adjust_counts) <- dimnames(counts)
  #return(adjust_counts)
  
  ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
  adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[keep, ] <- adjust_counts
  adjust_counts_whole[rm, ] <- countsOri[rm, ]
  return(adjust_counts_whole)
}

################################################################################################
################################################################################################
################################################################################################

# Example usage of ComBat_seq_custom


# # Import unnormalized raw counts as data frames
# # Extract flight and ground control samples, combine
# raw_242 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-242/GLDS-242_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# raw_242_flt <- as.matrix(raw_242[,grep("FLT", names(raw_242))])
# raw_242_gc <- as.matrix(raw_242[,grepl("GC", names(raw_242))])
# raw_242_flt_gc <- cbind(raw_242_flt, raw_242_gc)
# 
# raw_245 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-245/GLDS-245_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# raw_245_flt <- as.matrix(raw_245[,grepl("FLT", names(raw_245))])
# raw_245_gc <- as.matrix(raw_245[,grepl("GC", names(raw_245))])
# raw_245_flt_gc <- cbind(raw_245_flt, raw_245_gc)
# 
# raw_counts <- cbind(raw_242_flt_gc, raw_245_flt_gc)
# 
# bc_batchqc <- matrix(c(rep.int(242, 9),rep.int(245, 39)))
# # 1 = flight, 2 = ground control
# cov_batchqc <- matrix(c(rep.int(1, 5),rep.int(2, 4),rep.int(1, 20),rep.int(2, 19)))
# # batchQC(raw_counts, batch=bc_batchqc, report_file = "batchqc_flt_gc_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# # batchQC(raw_counts, batch=bc_batchqc, condition=cov_batchqc, report_file = "batchqc_flt_gc_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")



# bc_combat <- c(rep.int(242, 9),rep.int(245, 39))
# cov_combat <- c(rep.int(1, 5),rep.int(2, 4),rep.int(1, 20),rep.int(2, 19))
# # combat_seq_auto_counts <- ComBat_seq(counts=raw_counts, batch=bc_combat, group=cov_combat)


## Generate PCA for UMRR data
# 
# size_var <- 3
# alpha_var <- 0.6
# 
# pca_umrr_groups <- matrix(c(rep("242", 9), rep("245", 6)))
# rownames(pca_umrr_groups) <- colnames(raw_umrr)
# colnames(pca_umrr_groups) <- "Dataset"
# 
# combat_seq_umrr_in <- log2(ceiling(combat_seq_umrr))
# combat_seq_umrr_in[is.infinite(combat_seq_umrr_in)] <- 0
# combat_seq_umrr_in[is.nan(combat_seq_umrr_in)] <- 0
# 
# PCA_combat_seq_umrr <- prcomp(t(combat_seq_umrr_in), scale = FALSE)
# 
# autoplot(PCA_combat_seq_umrr, data=pca_umrr_groups, colour='Dataset', size=size_var, alpha=alpha_var) +
#   theme_classic()


## Run ComBat_seq
# combat_seq_man_counts <- ComBat_seq_custom(counts=raw_counts, batch=bc_combat, custom.gamma=umrr_gamma, custom.phi=umrr_phi, group=cov_combat)
# combat_seq_both_counts <- ComBat_seq(counts=combat_seq_man_counts, batch=bc_combat, group=cov_combat)


## Set all negative values to 0 for DESeq2 use
# raw_counts_dge <- raw_counts
# combat_seq_man_counts_dge <- combat_seq_man_counts
# combat_seq_auto_counts_dge <- combat_seq_auto_counts
# combat_seq_both_counts_dge <- combat_seq_both_counts
# raw_counts_dge[raw_counts_dge<0] <- 0
# combat_seq_man_counts_dge[combat_seq_man_counts_dge<0] <- 0
# combat_seq_auto_counts_dge[combat_seq_auto_counts_dge<0] <- 0
# combat_seq_both_counts_dge[combat_seq_both_counts_dge<0] <- 0

## Export counts
# Name: GLDS-242-245_FLT-GC
## write.csv(raw_counts_dge,file.path("~/Documents/SLSTP/Raw/GLDS-242-245_FLT-GC/Uncorrected", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(combat_seq_man_counts_dge,file.path("~/Documents/SLSTP/Raw/GLDS-242-245_FLT-GC/ComBat_seq-URR", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(combat_seq_auto_counts_dge,file.path("~/Documents/SLSTP/Raw/GLDS-242-245_FLT-GC/ComBat_seq-Standard", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(combat_seq_both_counts_dge,file.path("~/Documents/SLSTP/Raw/GLDS-242-245_FLT-GC/ComBat_seq-Both", "Unnormalized_Counts.csv"), row.names = TRUE)

## BatchQC tests
# batchQC(combat_seq_man_counts, batch=bc_batchqc, condition=cov_batchqc, report_file = "batchqc_flt_gc_combat_seq_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# batchQC(combat_man_counts, batch=bc_batchqc, condition=cov_batchqc, report_file = "batchqc_flt_gc_man_combat_seq_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# combat_seq_both_counts <- ComBat_seq(dat=combat_man_counts, batch=bc_combat, mod=cov_combat, par.prior=TRUE, prior.plots=FALSE)
# batchQC(combat_both_counts, batch=bc_batchqc, condition=cov_batchqc, report_file = "batchqc_flt_gc_both_combat_seq_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")

################################################################################################
################################################################################################
################################################################################################
