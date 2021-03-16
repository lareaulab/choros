load_bam <- function(bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                     f5_length=2, f3_length=3, full=F, nt_base=F, num_cores=NULL) {
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
  num_footprints <- nrow(alignment)
  print(paste("Read in", num_footprints, "total footprints"))
  print(paste("... Removing",
              sum(is.na(alignment$rname)),
              paste0("(", round(sum(is.na(alignment$rname)) / num_footprints * 100, 1), "%)"),
              "unaligned footprints"))
  alignment <- subset(alignment, !is.na(alignment$rname))
  # 2. assign 5' UTR lengths
  transcript_length <- load_lengths(transcript_length_fname)
  alignment$utr5_length <- transcript_length$utr5_length[match(alignment$rname,
                                                               transcript_length$transcript)]
  # 3. calculate frame
  alignment$frame <- (alignment$pos - alignment$utr5_length - 1) %% 3
  # 4. calculate 5' and 3' digest lengths
  offsets <- load_offsets(offsets_fname)
  alignment$d5 <- offsets$offset[prodlim::row.match(alignment[,c("frame", "qwidth")],
                                                    offsets[c("frame", "length")])]
  print(paste("... Removing",
              sum(is.na(alignment$d5)),
              paste0("(", round(sum(is.na(alignment$d5)) / num_footprints * 100, 1), "%)"),
              "footprints outside A site offset definitions"))
  alignment <- subset(alignment, !is.na(alignment$d5))
  # 5. calculate 3' digest lengths
  alignment$d3 <- with(alignment, qwidth - d5 - 3)
  # 6. calculate cod_idx, remove footprints mapped outside coding region
  alignment$cod_idx <- with(alignment, (pos + d5 - utr5_length + 2) / 3)
  alignment$cds_length <- transcript_length$cds_length[match(alignment$rname,
                                                             transcript_length$transcript)]/3
  outside_cds <- ((alignment$cod_idx <= 0) | (alignment$cod_idx > alignment$cds_length))
  print(paste("... Removing",
              sum(outside_cds),
              paste0("(", round(sum(outside_cds) / num_footprints * 100, 1), "%)"),
              "footprints outside CDS"))
  alignment <- subset(alignment, !outside_cds)
  # 7. pull bias sequences
  alignment$rpf_f5 <- substr(alignment$seq, 1, f5_length)
  alignment$rpf_f3 <- mapply(substr, alignment$seq, alignment$qwidth-f3_length+1, alignment$qwidth)
  invalid_bias_seq <- (grepl("N", alignment$rpf_f5) | grepl("N", alignment$rpf_f3))
  print(paste("... Removing",
              sum(invalid_bias_seq),
              paste0("(", round(sum(invalid_bias_seq) / num_footprints * 100, 1), "%)"),
              "footprints with N in bias region"))
  alignment <- subset(alignment, !invalid_bias_seq)
  alignment$rpf_f5 <- factor(alignment$rpf_f5)
  alignment$rpf_f3 <- factor(alignment$rpf_f3)
  # 8. aggregate alignments
  alignment <- aggregate(tag.ZW ~ rname + utr5_length + cod_idx + d5 + d3 + rpf_f5 + rpf_f3,
                         data=alignment, FUN=sum)
  # 8. return nt_base
  if(nt_base) {
    alignment$nt_base <- substr(as.character(alignment$rpf_f5), 1, 1)
    alignment$nt_base[!grepl("^0(A|T|C|G)", as.character(alignment$tag.MD))] <- "-"
    nt_bases <- c("-", "A", "T", "C", "G")
    alignment$nt_base <- factor(alignment$nt_base, levels=nt_bases)
    num_mismatch <- sum(alignment$nt_base != "-")
    print(paste("...",
                num_mismatch,
                paste0("(", round(num_mismatch / nrow(alignment) * 100, 1), "%)"),
                "footprints with non-templated 5' base"))
    for(nt in nt_bases[-1]) {
      print(paste("...",
                  sum(alignment$nt_base == nt),
                  paste0("(", round(sum(alignment$nt_base==nt) / num_mismatch * 100, 1), "%)")))
    }
    alignment$nt_d5 <- alignment$d5
    which_ntBase <- which(as.character(alignment$nt_base) != "-")
    alignment$nt_d5[which_ntBase] <- alignment$nt_d5[which_ntBase] - 1
  }
  # 9. return genome-annotated bias sequences
  chunks <- cut(seq.int(nrow(alignment)), num_cores)
  transcript_seq <- load_fa(transcript_fa_fname)
  alignment <- foreach(x=split(alignment, chunks),
                       .combine='rbind', .export=c("get_bias_seq")) %dopar% {
                         within(x, {
                           genome_f5 <- mapply(get_bias_seq,
                                               as.character(rname), cod_idx, d5, utr5_length,
                                               MoreArgs=list(transcript_seq=transcript_seq,
                                                             bias_region="f5",
                                                             bias_length=f5_length))
                           genome_f3 <- mapply(get_bias_seq,
                                               as.character(rname), cod_idx, d3, utr5_length,
                                               MoreArgs=list(transcript_seq=transcript_seq,
                                                             bias_region="f3",
                                                             bias_length=f3_length))
                         })
                       }
  # return data
  colnames(alignment)[colnames(alignment)=="rname"] <- "transcript"
  colnames(alignment)[colnames(alignment)=="tag.ZW"] <- "count"
  subset_features <- c("transcript", "cod_idx", "d5", "d3", "rpf_f5", "rpf_f3",
                       "genome_f5", "genome_f3", "count")
  if(nt_base) { subset_features <- c(subset_features, "nt_base", "nt_d5") }
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
  transcript_seq <- load_fa(transcript_fa_fname)
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
                                 row.names=NULL)
                    }
  if(!is.null(d5_d3_subsets)) {
    dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons, utr5_length),
                                   d5_d3_subsets)
  } else {
    dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons, utr5_length),
                                   expand.grid(d5=digest5_lengths, d3=digest3_lengths))
  }
  chunks <- cut(seq.int(nrow(dat)), num_cores)
  dat <- foreach(x=split(dat, chunks),
                 .combine='rbind', .export=c("get_bias_seq")) %dopar% {
                   within(x, {
                     genome_f5 <- mapply(get_bias_seq,
                                         as.character(transcript), cod_idx, d5, utr5_length,
                                         MoreArgs=list(transcript_seq=transcript_seq,
                                                       bias_region="f5",
                                                       bias_length=f5_length))
                     genome_f3 <- mapply(get_bias_seq,
                                         as.character(transcript), cod_idx, d3, utr5_length,
                                         MoreArgs=list(transcript_seq=transcript_seq,
                                                       bias_region="f3",
                                                       bias_length=f3_length))
                   })
                 }
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
  match_rows <- prodlim::row.match(regression_data[, features], bam_dat[, features])
  counts <- bam_dat[match_rows, which_column]
  counts[is.na(counts)] <- 0
  counts <- round(counts, digits=0) # return integer counts for glm.nb()
  return(counts)
}
