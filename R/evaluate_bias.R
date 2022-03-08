evaluate_bias <- function(dat, which_column="count",
                          transcripts_fa_fname, transcripts_length_fname,
                          trunc5=20, trunc3=20, num_f5_codons=6, num_f3_codons=6,
                          type="codon", metric="corr") {
  # perform iXnos regression and generate leave-one-out correlation plots
  ## dat: data.frame containing regression predictors and column of corrected counts
  ## which_column: character; name of column containing counts
  ## transcripts_fa_fname: character; filepath to transcriptome .fa file
  ## transcripts_length_fname: character; filepath to transcriptome lengths file
  ## trunc5: integer; number of 5' codons to leave out in codon correlation regression
  ## trunc3: integer; number of 3' codons to leave out in codon correlation regression
  ## num_f5_codons: integer; number of codons to include 5' of A site in regression
  ## num_f3_codons: integer; number of codons to include 3' of A site in regression
  ## type: character; one of "codon" or "nt"
  ## metric: character; one of:
  ##### corr: delta correlation between full model and leave-one-out model
  ##### norm: norm of regression coefficients at that position
  transcripts_seq <- read_fasta_as_codons(transcripts_fa_fname, transcripts_length_fname)
  transcripts_lengths <- load_lengths(transcripts_length_fname)
  # 1. aggregate counts by codon
  cts_by_codon <- aggregate(formula(paste(which_column, "~ transcript + cod_idx")),
                            data=dat, FUN=sum, na.rm=T)
  # 2. remove truncated codons, scale footprint counts by transcript mean
  transcripts <- levels(dat$transcript)
  cts_by_codon_trunc <- lapply(transcripts,
                               function(x) {
                                 tmp_subset <- subset(cts_by_codon, transcript==x)
                                 # remove trunc5 and trunc3 codons
                                 num_codons <- transcripts_lengths$cds_length[transcripts_lengths$transcript==x]/3
                                 tmp_subset <- subset(tmp_subset,
                                                      cod_idx > trunc5 & cod_idx <= (num_codons-trunc3))
                                 # scale per-codon counts by transcript mean
                                 tmp_subset[, which_column] <- tmp_subset[, which_column] / mean(tmp_subset[, which_column])
                                 # return truncated & scaled counts
                                 return(tmp_subset)
                               })
  cts_by_codon <- do.call(rbind, cts_by_codon_trunc)
  ##### remove transcripts with NA counts (mean = 0?)
  ## TODO: error message
  bad_transcript <- which(sapply(cts_by_codon_trunc, function(x) any(is.na(x[, which_column]))))
  cts_by_codon_trunc <- sapply(bad_transcript, function(x) cts_by_codon_trunc[[x]]$transcript)
  cts_by_codon <- subset(cts_by_codon, !(transcript %in% as.character(cts_by_codon_trunc)))
  # 3. make data.frame for iXnos regression
  codons <- data.frame(t(sapply(1:nrow(cts_by_codon),
                                function(x) {
                                  transcript <- as.character(cts_by_codon$transcript[x])
                                  cod_idx <- cts_by_codon$cod_idx[x] - 1 # transcripts_seq uses 0-based indexing
                                  Asite_index <- which(names(transcripts_seq[[transcript]])==as.character(cod_idx))
                                  codons <- transcripts_seq[[transcript]][(Asite_index-num_f5_codons):(Asite_index+num_f3_codons)]
                                  return(codons)
                                })), stringsAsFactors=F)
  colnames(codons) <- c(paste0("n", num_f5_codons:3),
                        "E", "P", "A",
                        paste0("p", 1:num_f3_codons))
  if(type=="nt") {
    codons <- data.frame(lapply(codons, as.character), stringsAsFactors=F)
    codons <- sapply(seq(nrow(codons)), function(x) {paste0(codons[x,], collapse="")})
    codons <- data.frame(matrix(unlist(strsplit(codons, split="")),
                                ncol=3*(num_f5_codons+num_f3_codons+1), byrow=T),
                         stringsAsFactors=F)
    colnames(codons) <- c(paste0("n", (3*num_f3_codons):7),
                          paste0(rep(c("E", "P", "A"), each=3), 0:2),
                          paste0("p", 1:(3*num_f5_codons)))
  }
  count_dat <- data.frame(count = cts_by_codon[, which_column], codons,
                          stringsAsFactors=F)
  if(metric=="corr") {
    # 4. full model
    count_full_cor <- cor(count_dat$count, predict(lm(count ~ ., data=count_dat)))
    # 5. leave-one-out models
    count_loo_cor <- sapply(colnames(codons),
                            function(x) {
                              cor(count_dat$count,
                                  predict(lm(count ~ ., data=count_dat[,-which(colnames(count_dat)==x)])))
                            })
    # 6. return model correlations
    # model_cor <- c(count_full_cor, count_loo_cor)
    # names(model_cor) <- c("full", colnames(codons))
    model_cor <- count_full_cor - count_loo_cor
    return(model_cor)
  } else {
    # 4. iXnos regression
    fit <- lm(count ~ ., data=count_dat)
    # 5. return norm of model coefficients per position
    model_norm <- sapply(colnames(codons),
                         function(x) {
                           sum(coef(fit)[grepl(x, names(coef(fit)))]^2)
                         })
    # model_norm <- model_norm / model_norm["A"]
    return(model_norm)
  }
}
