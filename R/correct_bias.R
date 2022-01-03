correct_bias_sim <- function(dat, p5bias, n3bias, which_column="count",
                             which_f5="genome_f5", which_f3="genome_f3") {
  # predict counts if bias sequences were set to reference level, using simulation parameters
  ## dat: data.frame containing regression predictors
  ## p5bias: named numeric vector; recovery probabilities from simulation
  ## n3bias: named numeric vector; recovery probabilities from simulation
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  # 1. establish f5 scaling factors
  f5_ref_level <- levels(dat[, which_f5])[1]
  f5_coefs <- p5bias / p5bias[f5_ref_level]
  f5_coefs[f5_coefs==0] <- 1
  # 2. establish f3 scaling factors
  f3_ref_level <- levels(dat[, which_f3])[1]
  f3_coefs <- n3bias / n3bias[f3_ref_level]
  f3_coefs[f3_coefs==0] <- 1
  # 3. calculate corrected counts
  corrected_count <- dat[, which_column] / (f5_coefs[as.character(dat[, which_f5])] *
                                              f3_coefs[as.character(dat[, which_f3])])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return dat with corrected counts
  return(corrected_count)
}

correct_bias_noInteraction <- function(dat, fit, which_column="count",
                                       which_f5="genome_f5", which_f3="genome_f3",
                                       fit_f5="genome_f5", fit_f3="genome_f3") {
  # predict counts if bias sequences were set to reference level, allowing residuals
  ## dat: data.frame containing regression predictors
  ## fit: glm object ; output from MASS::glm.nb()
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  ## fit_f5: character; prefix to f5 coefficients in fit object
  ## fit_f3: character; frefixprefix to f3 coefficients in fit object
  # 1. establish f5 scaling factors
  f5_coefs <- coef(fit)[match(paste0(fit_f5, fit$xlevels[[fit_f5]]), names(coef(fit)))]
  names(f5_coefs) <- fit$xlevels[[fit_f5]]
  f5_coefs[is.na(f5_coefs)] <- 0
  f5_coefs <- exp(f5_coefs)
  # 2. establish f3 scaling factors
  f3_coefs <- coef(fit)[match(paste0(fit_f3, fit$xlevels[[fit_f3]]), names(coef(fit)))]
  names(f3_coefs) <- fit$xlevels[[fit_f3]]
  f3_coefs[is.na(f3_coefs)] <- 0
  f3_coefs <- exp(f3_coefs)
  # 3. calculate corrected counts
  corrected_count <- dat[, which_column] / (f5_coefs[as.character(dat[, which_f5])] *
                                              f3_coefs[as.character(dat[, which_f3])])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return corrected counts
  return(corrected_count)
}

correct_bias <- function(dat, intrxn_fit, which_column="count",
                         which_f5="genome_f5", which_f3="genome_f3",
                         fit_f5="genome_f5", fit_f3="genome_f3") {
  # correct biased counts with coefficients from interaction regression model
  ## dat: data.frame containing regression predictors
  ## intrxn_fit: glm object ; output from MASS::glm.nb()
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  ## fit_f5: character; prefix to f5 coefficients in fit object
  ## fit_f3: character; frefixprefix to f3 coefficients in fit object
  fit_coefs <- coef(intrxn_fit)
  fit_coefs[abs(fit_coefs) > 10] <- 0
  # 1. establish f5 scaling factors
  f5_ref <- intrxn_fit$xlevels[[fit_f5]][1]
  d5_ref <- intrxn_fit$xlevels$d5[1]
  f5_coefs <- expand.grid(f5=intrxn_fit$xlevels[[fit_f5]], d5=intrxn_fit$xlevels$d5, stringsAsFactors=F)
  f5_coefs$correction <- sapply(seq(nrow(f5_coefs)),
                                function(x) {
                                  tmp_f5 <- f5_coefs$f5[x]
                                  tmp_d5 <- f5_coefs$d5[x]
                                  if(tmp_f5 == f5_ref) {
                                    tmp_coef <- 0
                                  } else {
                                    tmp_coef <- fit_coefs[paste0(fit_f5, tmp_f5)]
                                  }
                                  if(tmp_d5 != d5_ref & tmp_f5 != f5_ref) {
                                    tmp_intrxn <- paste0("d5", tmp_d5, ":", fit_f5, tmp_f5)
                                    if(!is.na(fit_coefs[tmp_intrxn])) {
                                      tmp_coef <- tmp_coef + fit_coefs[tmp_intrxn]
                                    }
                                  }
                                  return(tmp_coef)
                                })
  f5_coefs$correction <- exp(f5_coefs$correction)
  # 2. establish f3 scaling factors
  f3_ref <- intrxn_fit$xlevels[[fit_f3]][1]
  d3_ref <- intrxn_fit$xlevels$d3[1]
  f3_coefs <- expand.grid(f3=intrxn_fit$xlevels[[fit_f3]], d3=intrxn_fit$xlevels$d3, stringsAsFactors=F)
  f3_coefs$correction <- sapply(seq(nrow(f3_coefs)),
                                function(x) {
                                  tmp_f3 <- f3_coefs$f3[x]
                                  tmp_d3 <- f3_coefs$d3[x]
                                  if(tmp_f3 == f3_ref) {
                                    tmp_coef <- 0
                                  } else {
                                    tmp_coef <- fit_coefs[paste0(fit_f3, tmp_f3)]
                                  }
                                  if(tmp_d3 != d3_ref & tmp_f3 != f3_ref) {
                                    tmp_intrxn <- paste0("d3", tmp_d3, ":", fit_f3, tmp_f3)
                                    if(!is.na(fit_coefs[tmp_intrxn])) {
                                      tmp_coef <- tmp_coef + fit_coefs[tmp_intrxn]
                                    }
                                  }
                                  return(tmp_coef)
                                })
  f3_coefs$correction <- exp(f3_coefs$correction)
  # 3. calculate corrected counts
  f5_indices <- prodlim::row.match(dat[, c(which_f5, "d5")], f5_coefs[, c("f5", "d5")])
  f3_indices <- prodlim::row.match(dat[, c(which_f3, "d3")], f3_coefs[, c("f3", "d3")])
  corrected_count <- dat[, which_column] / (f5_coefs$correction[f5_indices] * f3_coefs$correction[f3_indices])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return corrected counts
  return(corrected_count)
}

correct_bias_frameBySize <- function(dat, nb_fits, which_column="count",
                                     which_f5="genome_f5", which_f3="genome_f3",
                                     fit_f5="genome_f5", fit_f3="genome_f3") {
  # correct counts for individual footprints
  ## dat: data.frame containing regression predictors
  ## nb_fits: list of glm objects ; names correspond to d5/d3 subsets
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  ## fit_f5: character; prefix to f5 coefficients in fit object
  ## fit_f3: character; prefix to f3 coefficients in fit object
  # 1. establish f5 scaling factors
  f5_coefs <- do.call(rbind,
                      lapply(seq_along(nb_fits),
                             function(x) {
                               subset_coefs <- coef(nb_fits[[x]])
                               subset_d5 <- as.numeric(strsplit(names(nb_fits)[x], split="_")[[1]][2])
                               subset_d3 <- as.numeric(strsplit(names(nb_fits)[x], split="_")[[1]][4])
                               tmp_coefs <- data.frame(d5=subset_d5,
                                                       d3=subset_d3,
                                                       f5=nb_fits[[x]]$xlevels[[fit_f5]],
                                                       stringsAsFactors=F)
                               tmp_coefs$correction <- subset_coefs[match(paste0(fit_f5, tmp_coefs$f5),
                                                                          names(subset_coefs))]
                               tmp_coefs$correction[is.na(tmp_coefs$correction)] <- 0
                               tmp_coefs$correction <- exp(tmp_coefs$correction)
                               return(tmp_coefs)
                             }))
  # 2. establish f3 scaling factors
  f3_coefs <- do.call(rbind,
                      lapply(seq_along(nb_fits),
                             function(x) {
                               subset_coefs <- coef(nb_fits[[x]])
                               subset_d5 <- as.numeric(strsplit(names(nb_fits)[x], split="_")[[1]][2])
                               subset_d3 <- as.numeric(strsplit(names(nb_fits)[x], split="_")[[1]][4])
                               tmp_coefs <- data.frame(d5=subset_d5,
                                                       d3=subset_d3,
                                                       f3=nb_fits[[x]]$xlevels[[fit_f3]],
                                                       stringsAsFactors=F)
                               tmp_coefs$correction <- subset_coefs[match(paste0(fit_f3, tmp_coefs$f3),
                                                                          names(subset_coefs))]
                               tmp_coefs$correction[is.na(tmp_coefs$correction)] <- 0
                               tmp_coefs$correction <- exp(tmp_coefs$correction)
                               return(tmp_coefs)
                             }))
  # 3. calculate corrected counts
  f5_indices <- prodlim::row.match(dat[, c("d5", "d3", which_f5)], f5_coefs[, c("d5", "d3", "f5")])
  f3_indices <- prodlim::row.match(dat[, c("d5", "d3", which_f3)], f3_coefs[, c("d5", "d3", "f3")])
  corrected_count <- dat[, which_column] / (f5_coefs$correction[f5_indices] * f3_coefs$correction[f3_indices])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return corrected counts
  return(corrected_count)
}

correct_bias_zinb <- function(dat, fit, which_column="count",
                              which_f5="genome_f5", which_f3="genome_f3",
                              fit_f5="mod_f5", fit_f3="genome_f3") {
  # predict counts if bias sequences were set to reference level, allowing residuals
  ## dat: data.frame containing regression predictors
  ## fit: glm object ; output from pscl::zeroinfl()
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  ## fit_f5: character; prefix to f5 coefficients in fit object
  ## fit_f3: character; frefixprefix to f3 coefficients in fit object
  # 1. establish f5 scaling factors
  f5_coefs <- coef(fit)[match(paste0("count_f5", levels(fit$model[[fit_f5]])), names(coef(fit)))]
  names(f5_coefs) <- levels(fit$model[[fit_f5]])
  f5_coefs[is.na(f5_coefs)] <- 0
  f5_coefs <- exp(f5_coefs)
  # 2. establish f3 scaling factors
  f3_coefs <- coef(fit)[match(paste0("count_f3", levels(fit$model[[fit_f3]])), names(coef(fit)))]
  names(f3_coefs) <- levels(fit$model[[fit_f3]])
  f3_coefs[is.na(f3_coefs)] <- 0
  f3_coefs <- exp(f3_coefs)
  # 3. calculate corrected counts
  corrected_count <- dat[, which_column] / (f5_coefs[as.character(dat[, which_f5])] *
                                              f3_coefs[as.character(dat[, which_f3])])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return corrected counts
  return(corrected_count)
}

correct_bias_nt <- function(dat, fit, which_column="count", which_bias=NULL) {
  # predict counts if bias sequences were set to reference level, allowing residuals
  ## dat: data.frame containing regression predictors
  ## fit: glm object ; output from pscl::zeroinfl()
  ## which_column: character; name of column containing uncorrected counts
  ## which_bias: character vector; column names corresponding to end biases to be corrected
  # 1. establish scaling factors
  fit_coefs <- coef(fit)
  if(is.null(which_bias)) {
    which_bias <- grep("^f", colnames(dat), value=T)
    which_bias <- which_bias[which_bias %in% names(fit$xlevels)]
  }
  scaling_factors <- lapply(which_bias,
                            function(x) {
                              tmp_coefs <- fit_coefs[match(paste0(x, fit$xlevels[[x]]), names(fit_coefs))]
                              names(tmp_coefs) <- fit$xlevels[[x]]
                              tmp_coefs[is.na(tmp_coefs)] <- 0
                              tmp_coefs <- exp(tmp_coefs)
                            })
  names(scaling_factors) <- which_bias
  # 2. calculate corrected counts
  correction <- lapply(which_bias,
                       function(x) {
                         scaling_factors[[x]][as.character(dat[, x])]
                       })
  correction <- do.call(cbind, correction)
  correction <- exp(rowSums(log(correction)))
  corrected_count <- dat[, which_column] / correction
  # 3. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 4. return corrected counts
  return(corrected_count)
}

correct_bias_f5Interaction <- function(dat, intrxn_fit, which_column="count",
                                       which_f5="genome_f5", which_f3="genome_f3",
                                       fit_f5="genome_f5", fit_f3="genome_f3") {
  # correct biased counts with coefficients from interaction regression model
  ## dat: data.frame containing regression predictors
  ## intrxn_fit: glm object ; output from MASS::glm.nb()
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  ## fit_f5: character; prefix to f5 coefficients in fit object
  ## fit_f3: character; prefix to f3 coefficients in fit object
  fit_coefs <- coef(intrxn_fit)
  fit_coefs[abs(fit_coefs) > 10] <- 0
  # 1. establish f5 scaling factors
  f5_ref <- intrxn_fit$xlevels[[fit_f5]][1]
  d5_ref <- intrxn_fit$xlevels$d5[1]
  f5_coefs <- expand.grid(f5=intrxn_fit$xlevels[[fit_f5]], d5=intrxn_fit$xlevels$d5, stringsAsFactors=F)
  f5_coefs$correction <- sapply(seq(nrow(f5_coefs)),
                                function(x) {
                                  tmp_f5 <- f5_coefs$f5[x]
                                  tmp_d5 <- f5_coefs$d5[x]
                                  if(tmp_f5 == f5_ref) {
                                    tmp_coef <- 0
                                  } else {
                                    tmp_coef <- fit_coefs[paste0(fit_f5, tmp_f5)]
                                  }
                                  if(tmp_d5 != d5_ref & tmp_f5 != f5_ref) {
                                    tmp_intrxn <- paste0("d5", tmp_d5, ":", fit_f5, tmp_f5)
                                    if(!is.na(fit_coefs[tmp_intrxn])) {
                                      tmp_coef <- tmp_coef + fit_coefs[tmp_intrxn]
                                    }
                                  }
                                  return(tmp_coef)
                                })
  f5_coefs$correction <- exp(f5_coefs$correction)
  # 2. establish f3 scaling factors
  f3_coefs <- coef(intrxn_fit)[match(paste0(fit_f3, intrxn_fit$xlevels[[fit_f3]]), names(coef(intrxn_fit)))]
  names(f3_coefs) <- intrxn_fit$xlevels[[fit_f3]]
  f3_coefs[is.na(f3_coefs)] <- 0
  f3_coefs <- exp(f3_coefs)
  # 3. calculate corrected counts
  f5_indices <- prodlim::row.match(dat[, c(which_f5, "d5")], f5_coefs[, c("f5", "d5")])
  corrected_count <- dat[, which_column] / (f5_coefs$correction[f5_indices] * f3_coefs[as.character(dat[, which_f3])])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return corrected counts
  return(corrected_count)
}

correct_bias.glmnet <- function(dat, glmnet_fit, lambda=NULL, which_column="count",
                                which_f5="genome_f5", which_f3="genome_f3",
                                fit_f5="genome_f5", fit_f3="genome_f3") {
  # correct biased counts with coefficients from interaction regression model
  ## dat: data.frame containing regression predictors
  ## glmnet_fit: glmnet object ; output from glmnet())
  ## lambda: numeric; value of penalty parameter lambda to evaluate model
  ## which_column: character; name of column containing uncorrected counts
  ## which_f5: character; name of column containing 5' bias sequences to be used in correction
  ## which_f3: character; name of column containing 3' bias sequences to be used in correction
  ## fit_f5: character; prefix to f5 coefficients in fit object
  ## fit_f3: character; prefix to f3 coefficients in fit object
  if(is.null(lambda)) {
    lambda <- glmnet_fit$lambda[which.max(glmnet_fit$dev.ratio)]
  }
  fit_coefs <- coef(glmnet_fit, s=lambda)[,1]
  xlevels <- unlist(glmnet_fit$xlev, recursive=F)
  # 1. establish f5 scaling factors
  f5_ref <- xlevels[[fit_f5]][1]
  d5_ref <- xlevels$d5[1]
  f5_coefs <- expand.grid(f5=xlevels[[fit_f5]], d5=xlevels$d5, stringsAsFactors=F)
  f5_coefs$correction <- sapply(seq(nrow(f5_coefs)),
                                function(x) {
                                  tmp_f5 <- f5_coefs$f5[x]
                                  tmp_d5 <- f5_coefs$d5[x]
                                  if(tmp_f5 == f5_ref) {
                                    tmp_coef <- 0
                                  } else {
                                    tmp_coef <- fit_coefs[paste0(fit_f5, tmp_f5)]
                                  }
                                  if(tmp_d5 != d5_ref & tmp_f5 != f5_ref) {
                                    tmp_glmnet <- paste0("d5", tmp_d5, ":", fit_f5, tmp_f5)
                                    if(!is.na(fit_coefs[tmp_glmnet])) {
                                      tmp_coef <- tmp_coef + fit_coefs[tmp_glmnet]
                                    }
                                  }
                                  return(tmp_coef)
                                })
  f5_coefs$correction <- exp(f5_coefs$correction)
  # 2. establish f3 scaling factors
  f3_ref <- xlevels[[fit_f3]][1]
  d3_ref <- xlevels$d3[1]
  f3_coefs <- expand.grid(f3=xlevels[[fit_f3]], d3=xlevels$d3, stringsAsFactors=F)
  f3_coefs$correction <- sapply(seq(nrow(f3_coefs)),
                                function(x) {
                                  tmp_f3 <- f3_coefs$f3[x]
                                  tmp_d3 <- f3_coefs$d3[x]
                                  if(tmp_f3 == f3_ref) {
                                    tmp_coef <- 0
                                  } else {
                                    tmp_coef <- fit_coefs[paste0(fit_f3, tmp_f3)]
                                  }
                                  if(tmp_d3 != d3_ref & tmp_f3 != f3_ref) {
                                    tmp_glmnet <- paste0("d3", tmp_d3, ":", fit_f3, tmp_f3)
                                    if(!is.na(fit_coefs[tmp_glmnet])) {
                                      tmp_coef <- tmp_coef + fit_coefs[tmp_glmnet]
                                    }
                                  }
                                  return(tmp_coef)
                                })
  f3_coefs$correction <- exp(f3_coefs$correction)
  # 3. calculate corrected counts
  f5_indices <- prodlim::row.match(dat[, c(which_f5, "d5")], f5_coefs[, c("f5", "d5")])
  f3_indices <- prodlim::row.match(dat[, c(which_f3, "d3")], f3_coefs[, c("f3", "d3")])
  corrected_count <- dat[, which_column] / (f5_coefs$correction[f5_indices] * f3_coefs$correction[f3_indices])
  # 4. rescale predicted counts so they sum to original footprint count
  corrected_count <- corrected_count * sum(dat[, which_column]) / sum(corrected_count, na.rm=T)
  # 5. return corrected counts
  return(corrected_count)
}
