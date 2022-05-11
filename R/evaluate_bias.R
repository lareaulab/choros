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
  transcripts_lengths$num_codons <- with(transcripts_lengths, cds_length/3)
  # 1. aggregate counts by codon
  codon_cts <- aggregate(formula(paste(which_column, "~ transcript + cod_idx")),
                         data=dat, FUN=sum, na.rm=T)
  # remove truncated positions, flesh out unobserved positions
  transcripts <- levels(dat$transcript)
  cts_by_codon <- lapply(transcripts,
                         function(x) {
                           num_codons <- transcripts_lengths$num_codons[transcripts_lengths$transcript==x]
                           return(data.frame(transcript=x,
                                             cod_idx=seq(trunc5+1, num_codons-trunc3)))
                         })
  cts_by_codon <- do.call(rbind, cts_by_codon)
  cts_by_codon$count <- codon_cts[[which_column]][match_rows(cts_by_codon, codon_cts)]
  cts_by_codon$count[is.na(cts_by_codon$count)] <- 0
  # 2. scale footprint counts by transcript mean
  cts_by_codon <- split(cts_by_codon, cts_by_codon$transcript)
  cts_by_codon <- lapply(cts_by_codon,
                         function(x) {
                           if(sum(x$count) == 0) {
                             print("ERROR: no footprints observed!")
                             return(NULL)
                           } else {
                             within(x, count <- count/mean(count, na.rm=T))
                           }
                         })
  cts_by_codon <- do.call(rbind, cts_by_codon)
  # 3. make data.frame for iXnos regression
  codons <- data.frame(t(sapply(seq(nrow(cts_by_codon)),
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
  count_dat <- data.frame(count = cts_by_codon$count, codons,
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
