load_bam <- function(bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                     f5_length=3, f3_length=3, full=F, num_cores=NULL) {
  # calculate proportion of footprints within each 5' and 3' digest length combination
  ## bam_fname: character; file.path to .bam alignment file
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## f5_length: integer; length of 5' bias region
  ## f3_length: integer; length of 3' bias region
  ## full: logical; whether to import all fields in .bam alignment file
  ## num_cores: integer; number of cores to parallelize over
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)
  # 1. read in footprints
  bam_file <- Rsamtools::BamFile(bam_fname)
  if(full) {
    features <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual", "qwidth")
  } else {
    features <- c("rname", "pos", "seq", "qwidth")
  }
  bam_param <- Rsamtools::ScanBamParam(tag=c("ZW", "MD"), what=features)
  alignment <- data.frame(Rsamtools::scanBam(bam_file, param=bam_param)[[1]])
  num_footprints <- sum(alignment$tag.ZW, na.rm=T)
  print(paste("Read in", round(num_footprints, 1), "total RPF counts"))
  tmp_counts <- sum(subset(alignment, is.na(rname))$tag.ZW, na.rm=T)
  print(paste("... Removing",
              round(tmp_counts, 1),
              paste0("(", round(tmp_counts / num_footprints * 100, 1), "%)"),
              "unaligned RPF counts"))
  alignment <- subset(alignment, !is.na(alignment$rname))
  # 2. assign 5' UTR lengths
  transcript_length <- load_lengths(transcript_length_fname)
  alignment$utr5_length <- transcript_length$utr5_length[match(alignment$rname,
                                                               transcript_length$transcript)]
  # 3. calculate frame
  alignment$frame <- (alignment$pos - alignment$utr5_length - 1) %% 3
  # 4. calculate 5' and 3' digest lengths
  offsets <- load_offsets(offsets_fname)
  tmp_counts <- sum(subset(alignment, !(qwidth %in% unique(offsets$length)))$tag.ZW, na.rm=T)
  print(paste("... Removing",
              round(tmp_counts, 1),
              paste0("(", round(tmp_counts / num_footprints * 100, 1), "%)"),
              "RPF counts outside length distribution"))
  alignment <- subset(alignment, qwidth %in% unique(offsets$length))
  alignment$d5 <- offsets$offset[match_rows(alignment, offsets,
                                            c("frame", "qwidth"),
                                            c("frame", "length"))]
  tmp_counts <- sum(subset(alignment, is.na(d5))$tag.ZW)
  print(paste("... Removing",
              round(tmp_counts, 1),
              paste0("(", round(tmp_counts / num_footprints * 100, 1), "%)"),
              "footprints outside A site offset definitions"))
  alignment <- subset(alignment, !is.na(alignment$d5))
  # 5. calculate 3' digest lengths
  alignment$d3 <- with(alignment, qwidth - d5 - 3)
  # 6. calculate cod_idx, remove footprints mapped outside coding region
  alignment$cod_idx <- with(alignment, (pos + d5 - utr5_length + 2) / 3)
  alignment$cds_length <- transcript_length$cds_length[match(alignment$rname,
                                                             transcript_length$transcript)]/3
  outside_cds <- ((alignment$cod_idx <= 0) | (alignment$cod_idx > alignment$cds_length))
  tmp_counts <- sum(subset(alignment, outside_cds)$tag.ZW, na.rm=T)
  print(paste("... Removing",
              round(tmp_counts, 1),
              paste0("(", round(tmp_counts / num_footprints * 100, 1), "%)"),
              "footprints outside CDS"))
  alignment <- subset(alignment, !outside_cds)
  # 7. aggregate alignments
  alignment <- aggregate(tag.ZW ~ rname + utr5_length + cod_idx + d5 + d3,
                         data=alignment, FUN=sum)
  colnames(alignment)[colnames(alignment)=="rname"] <- "transcript"
  colnames(alignment)[colnames(alignment)=="tag.ZW"] <- "count"
  # 8. return genome-annotated bias sequences
  chunks <- cut(seq.int(nrow(alignment)), num_cores)
  transcript_seq <- load_fasta(transcript_fa_fname)
  alignment <- foreach(x=split(alignment, chunks),
                       .combine='rbind', .export=c("get_bias_seq")) %dopar% {
                         within(x, {
                           genome_f5 <- get_bias_seq(x, transcript_seq, "f5", f5_length)
                           genome_f3 <- get_bias_seq(x, transcript_seq, "f3", f3_length)
                         })
                       }
  # return data
  subset_features <- c("transcript", "cod_idx", "d5", "d3", "rpf_f5", "rpf_f3",
                       "genome_f5", "genome_f3", "count")
  if(!full) { alignment <- alignment[, subset_features] }
  return(alignment)
}

init_data <- function(transcript_fa_fname, transcript_length_fname,
                      digest5_lengths=15:18, digest3_lengths=9:11,
                      d5_d3_subsets=NULL, f5_length=2, f3_length=3,
                      num_cores=NULL, which_transcripts=NULL,
                      exclude_codons5=10, exclude_codons3=10) {
  # initialize data.frame for downstream GLM
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## digest5_lengths: integer vector; legal 5' digest lengths
  ## digest3_lengths: integer vector; legal 3' digest lengths
  ## d5_d3_subsets: data.frame; columns "d5" and "d3" of d5/d3 subsets to initiate data over
  ## bias_length: integer; length of bias sequence
  ## num_cores: integer; number of cores to parallelize over
  ## which_transcripts: character vector; transcripts selected for regression
  ## exclude_codons5: integer; number of codons to exclude from 5' end of transcript
  ## exclude_codons3: integer; number of codons to exclude from 3' end of transcript
  transcript_seq <- load_fasta(transcript_fa_fname)
  transcript_length <- load_lengths(transcript_length_fname)
  if(!is.null(which_transcripts)) {
    transcript_seq <- transcript_seq[which_transcripts]
    transcript_length <- subset(transcript_length, transcript %in% which_transcripts)
  }
  transcript <- unlist(mapply(rep, x=transcript_length$transcript,
                              times=(transcript_length$cds_length/3 -
                                       exclude_codons5 - exclude_codons3)))
  cod_idx <- unlist(lapply(transcript_length$cds_length/3,
                           function(x) {
                             seq(exclude_codons5 + 1, x - exclude_codons3)
                           }))
  utr5_length <- transcript_length$utr5_length[match(transcript,
                                                     transcript_length$transcript)]
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)
  codons <- data.frame(transcript=transcript,
                       cod_idx=cod_idx,
                       utr5_length=utr5_length,
                       stringsAsFactors=F)
  chunks <- cut(seq.int(nrow(codons)), num_cores)
  codons <- foreach(x=split(codons, chunks),
                    .combine='rbind', .export=c("get_codons")) %dopar% {
                      data.frame(t(with(x,
                                        mapply(get_codons,
                                               transcript, cod_idx, utr5_length,
                                               MoreArgs=list(transcript_seq=transcript_seq)))),
                                 row.names=NULL, stringsAsFactors=F)
                    }
  if(!is.null(d5_d3_subsets)) {
    dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons, utr5_length),
                                   d5_d3_subsets, stringsAsFactors=F)
  } else {
    dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons, utr5_length,
                                              stringsAsFactors=F),
                                   expand.grid(d5=digest5_lengths, d3=digest3_lengths))
  }
  chunks <- cut(seq.int(nrow(dat)), num_cores)
  dat <- foreach(x=split(dat, chunks),
                 .combine='rbind', .export=c("get_bias_seq")) %dopar% {
                   within(x, {
                     genome_f5 <- get_bias_seq(x, transcript_seq, "f5", f5_length)
                     genome_f3 <- get_bias_seq(x, transcript_seq, "f3", f3_length)
                   })
                 }
  dat$transcript <- as.factor(dat$transcript)
  dat$d5 <- factor(dat$d5, levels=unique(d5_d3_subsets$d5))
  dat$d3 <- factor(dat$d3, levels=unique(d5_d3_subsets$d3))
  dat$genome_f5 <- as.factor(dat$genome_f5)
  dat$genome_f3 <- as.factor(dat$genome_f3)
  dat$count <- 0
  return(dat)
}

count_d5_d3 <- function(bam_dat, plot_title="") {
  # count number of footprints per d5/d3 combination
  ## bam_dat: data.frame; output from load_bam()
  ## plot_title: character; plot title for output heatmap
  subset_count <- aggregate(count ~ d5 + d3, data=bam_dat, FUN=sum)
  subset_count <- subset_count[order(subset_count$count, decreasing=T),]
  subset_count$proportion <- sapply(seq(nrow(subset_count)),
                                    function(x) {
                                      sum(subset_count$count[1:x])/sum(subset_count$count)
                                    })
  subset_count_plot <- ggplot(subset_count, aes(x=d5, y=d3, fill=count)) + geom_tile(col="black") +
    scale_fill_gradient(low="white", high="blue", name="Count") + theme_classic() +
    geom_text(aes(label=paste0(round(count/sum(count)*100, 1), "%"))) +
    ggtitle(plot_title) + xlab("5' digestion length") + ylab("3' digestion length")
  return(list(counts=subset_count, plot=subset_count_plot))
}

count_footprints <- function(bam_dat, regression_data, which_column="count", nt_base=F) {
  # count up footprints by transcript, A site, and digest lengths
  ## bam_dat: data.frame; output from load_bam()
  ## regression_data: data.frame; output from init_data()
  ## which_column: character; name of column containing counts
  ## nt_base: logical; whether to account for non-templated bases
  # count up footprints
  bam_dat <- subset(bam_dat, transcript %in% levels(regression_data$transcript))
  if(nt_base) {
    features <- c("transcript", "cod_idx", "mod_d5", "d3", "nt_base")
    bam_dat <- aggregate(formula(paste(which_column, "~ transcript + cod_idx + mod_d5 + d3 + nt_base")),
                         data=bam_dat, FUN=sum, na.rm=T)
  } else {
    features <- c("transcript", "cod_idx", "d5", "d3")
    bam_dat <- aggregate(formula(paste(which_column, "~ transcript + cod_idx + d5 + d3")),
                         data=bam_dat, FUN=sum, na.rm=T)
  }
  # add counts to regression data.frame
  counts <- bam_dat[match_rows(regression_data, bam_dat, features), which_column]
  counts[is.na(counts)] <- 0
  counts <- round(counts, digits=0) # return integer counts for glm.nb()
  return(counts)
}

calculate_codon_density <- function(bam_dat, transcript_length,
                                    exclude_codons5=10, exclude_codons3=10,
                                    normalize=F, which_column="count") {
  # return vector of counts per codon (for individual transcript)
  ## bam_dat: data.frame; output (or subset of) from load_bam()
  ## transcript_length: integer; number of codons in transcript
  ## exclude_codons5: integer; number of codons to exclude from 5' end of transcript
  ## exclude_codons3: integer; number of codons to exclude from 3' end of transcript
  ## normalize: logical; whether to normalize codon counts by average across transcript
  # 1. initialize empty vector
  counts <- rep(0, transcript_length)
  names(counts) <- seq(transcript_length)
  # 2. count up footprints per codon
  bam_dat_cts <- aggregate(count ~ cod_idx, data=bam_dat, FUN=sum)
  # 3. populate counts vector
  counts[bam_dat_cts$cod_idx] <- bam_dat_cts$count
  # 4. remove first and last codons
  counts <- counts[(1+exclude_codons5):(length(counts)-exclude_codons3)]
  # 5. normalize by average across transcript
  if(normalize) {
    counts / mean(counts)
  }
  return(counts)
}

calculate_transcript_density <- function(bam_dat, transcript_length_fname,
                                         statistic=mean,
                                         exclude_codons5=10, exclude_codons3=10) {
  # compute mean/median footprint density per codon across transcript
  ## bam_dat: data.frame; output from load_bam()
  ## transcript_length_fname: character; file.path to transcriptome lengths file
  ## statistic: character; function name (ex. "mean", "median")
  ## exclude_codons5: integer; number of codons to exclude from 5' end of transcript
  ## exclude_codons3: integer; number of codons to exclude from 3' end of transcript
  transcript_lengths <- load_lengths(transcript_length_fname)
  # 1. subset by transcript
  bam_dat$transcript <- droplevels(bam_dat$transcript)
  per_transcript <- split(bam_dat, bam_dat$transcript)
  num_codons <- transcript_lengths$cds_length[match(names(per_transcript),
                                                    transcript_lengths$transcript)] / 3
  # 2. aggregate footprints by cod_idx
  per_transcript <- lapply(seq_along(per_transcript),
                           function(x) {
                             calculate_codon_density(per_transcript[[x]],
                                                     num_codons[x],
                                                     exclude_codons5, exclude_codons3)
                           })
  names(per_transcript) <- levels(bam_dat$transcript)
  # 3. calculate mean/median codon density per transcript
  per_transcript <- sapply(per_transcript, FUN=statistic)
  per_transcript <- sort(per_transcript, decreasing=T)
  return(per_transcript)
}
