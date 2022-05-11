correct_bias <- function(dat, f5_method="correct_interxn",
                         f3_method="correct_interxn", correct_gc=T,
                         which_column="count", fit_coefs) {
  # wrapper function to correct biases
  ## dat: data.frame containing regression predictors
  ## f5_method: character; name of function to apply for 5' correction
  ## f3_method: character; name of function to apply for 3' correction
  ## correct_gc: logical; whether to correct for RPF %GC
  ## which_column: character; name of column in dat containing uncorrected counts
  ## fit_coefs: data.frame; output from parse_coefs()
  rpf_count <- sum(dat[, which_column], na.rm=T)
  # 1. correct 5' bias
  ### TODO check for required arguments for f5_method
  dat$corrected <- do.call(f5_method,
                           args=list(dat=dat, which_column=which_column,
                                     which_region="f5", fit_coefs=fit_coefs))
  # 2. correct 3' bias
  ### TODO check for required arguments for f3_method
  dat$corrected <- do.call(f3_method,
                           args=list(dat=dat, which_column="corrected",
                                     which_region="f3", fit_coefs=fit_coefs))
  # 3. correct for %GC
  if(correct_gc) {
    ### TODO check that fit_coefs specified
    dat$corrected <- correct_gc(dat, fit_coefs, which_column="corrected")
  }
  # 4. rescale corrected counts so they sum to original footprint count
  dat$corrected <- dat$corrected * rpf_count / sum(dat$corrected, na.rm=T)
  # 5. return corrected counts
  return(dat$corrected)
}

correct_fixed <- function(dat, which_column, which_region, bias_params) {
  # use simulation parameters as scaling factors for bias correction
  ## dat: data.frame containing regression predictors
  ## which_column: character; name of column in dat containing uncorrected counts
  ## which_region: character; one of "f5" or "f3" corresponding to RPF region
  ## bias_params: named numeric vector; recovery probabilities from simulation
  # 1. establish scaling factors
  ref_level <- levels(dat[, which_region])[1]
  scaling_factors <- bias_params / bias_params[ref_level]
  scaling_factors[scaling_factors==0] <- 1
  # 2. calculate corrected counts
  corrected_count <- dat[, which_column] /
    scaling_factors[as.character(dat[, which_region])]
  # 3. return corrected counts
  return(corrected_count)
}

correct_marginal <- function(dat, which_column, which_region, fit_coefs) {
  # use regression coefficients as scaling factors for bias correction
  ## only use marginal effects (i.e. not interaction terms)
  ## dat: data.frame containing regression predictors
  ## which_column: character; name of column in dat containing uncorrected counts
  ## which_region: character; one of "f5" or "f3" corresponding to RPF region
  ## fit_coefs: data.frame; output from parse_coefs()
  # 1. establish scaling factors
  scaling_factors <- subset(fit_coefs, group==which_region)$estimate
  names(scaling_factors) <- subset(fit_coefs, group==which_region)$term
  scaling_factors <- exp(scaling_factors)
  # 2. calculate corrected counts
  corrected_count <- dat[, which_column] /
    scaling_factors[as.character(dat[, which_region])]
  # 3. return corrected counts
  return(corrected_count)
}

correct_interxn <- function(dat, which_column, which_region, fit_coefs) {
  # use regression coefficients as scaling factors for bias correction
  ## use marginal and interaction effects
  ## dat: data.frame containing regression predictors
  ## which_column: character; name of column in dat containing uncorrected counts
  ## which_region: character; one of "f5" or "f3" corresponding to RPF region
  ## fit_coefs: data.frame; output from parse_coefs()
  which_digest <- sub("f", "d", which_region)
  which_interxn <- paste0(which_digest, ":", which_region)
  # 1. establish scaling factors
  bias_terms <- subset(fit_coefs, group==which_region)$term
  digest_terms <- subset(fit_coefs, group==which_digest)$term
  scaling_factors <- expand.grid(bias_term=bias_terms, digest_term=digest_terms,
                                 stringsAsFactors=F)
  scaling_factors$correction <- sapply(seq(nrow(scaling_factors)),
                                       function(x) {
                                         tmp_bias <- scaling_factors$bias_term[x]
                                         tmp_digest <- scaling_factors$digest_term[x]
                                         tmp_digest <- paste0(which_digest, "=", tmp_digest)
                                         marginal_effect <- subset(fit_coefs,
                                                                   group==which_region &
                                                                     term==tmp_bias)$estimate
                                         interxn_effect <- subset(fit_coefs,
                                                                  group==which_interxn &
                                                                    term==tmp_bias &
                                                                    group_2==tmp_digest)$estimate
                                         if(length(interxn_effect) == 0) {
                                           return(marginal_effect)
                                         } else {
                                           return(marginal_effect + interxn_effect)
                                         }
                                       })
  scaling_factors$correction <- exp(scaling_factors$correction)
  # 2. calculate corrected counts
  dat_indices <- match_rows(dat, scaling_factors,
                            c(which_region, which_digest),
                            c("bias_term", "digest_term"))
  corrected_count <- dat[, which_column] /
    scaling_factors$correction[dat_indices]
  # 3. return corrected counts
  return(corrected_count)
}

correct_gc <- function(dat, which_column, fit_coefs) {
  # use regression coefficients as scaling factor for %GC correction
  ## dat: data.frame containing regression predictors
  ## which_column: character; name of column in dat containing uncorrected counts
  ## fit_coefs: data.frame; output from parse_coefs()
  # 1. establish scaling factors
  scaling_factor <- subset(fit_coefs, group=="gc")$estimate
  # 2. calculate corrected counts
  mean_gc <- sum(dat[, which_column]*dat$gc, na.rm=T) / sum(dat[, which_column], na.rm=T)
  corrected_count <- exp(log(dat[, which_column]) + scaling_factor*mean_gc - scaling_factor*dat$gc)
  # 3. return corrected counts
  return(corrected_count)
}
