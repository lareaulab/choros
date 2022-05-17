load_lengths <- function(lengths_fname) {
  # load transcript lengths table
  ## length_fname: character; file path to transcript lengths file
  transcript_lengths <- read.table(lengths_fname, stringsAsFactors=F,
                                   col.names=c("transcript", "utr5_length", "cds_length", "utr3_length"))
  return(transcript_lengths)
}

gff_to_lengths <- function(gff_fname) {
  # load gff annotation file of gene regions
  ## gff_fname: character; file path to gff annotation file
  ## gff file requirements:
  ## - first column is transcript name
  ## - third column is one of UTR5, CDS, UTR3
  gff <- read.table(gff_fname, col.names=c("seqid", "source", "type", "start", "end",
                                           "score", "strand", "phase", "attributes"),
                    stringsAsFactors=F)
  transcripts <- unique(gff$seqid)
  transcript_lengths <- data.frame(transcript = transcripts,
                                   utr5_length = sapply(transcripts,
                                                        function(x) {
                                                          dat <- subset(gff, seqid==x & type=="UTR5")
                                                          return(dat$end - dat$start + 1)
                                                        }),
                                   cds_length = sapply(transcripts,
                                                       function(x) {
                                                         dat <- subset(gff, seqid==x & type=="CDS")
                                                         return(dat$end - dat$start + 1)
                                                       }),
                                   utr3_length = sapply(transcripts,
                                                        function(x) {
                                                          dat <- subset(gff, seqid==x & type=="UTR3")
                                                          return(dat$end - dat$start + 1)
                                                        }),
                                   row.names=NULL, stringsAsFactors=F)
  return(transcript_lengths)
}

load_fasta <- function(transcript_fa_fname) {
  # load transcript sequences from genome .fa file
  ## transcripts_fa_fname: character; file path to transcriptome .fa file
  transcript_sequences <- Biostrings::readDNAStringSet(transcript_fa_fname)
  transcript_sequences <- as.character(transcript_sequences)
  names(transcript_sequences) <- sapply(names(transcript_sequences),
                                        function(x) {
                                          strsplit(x, split=" ")[[1]][1]
                                        })
  return(transcript_sequences)
}

write_fasta <- function(transcript_seq, fa_fname) {
  # write transcript sequence to .fa file
  ## transcript_seq: character vector; transcript sequences
  ### - names of transcript_seq will become header for transcript sequences
  ## fa_fname: character; filepath to output .fa file
  transcript_seq <- Biostrings::DNAStringSet(transcript_seq)
  Biostrings::writeXStringSet(transcript_seq, fa_fname)
}

load_offsets <- function(offsets_fname) {
  # load A site offset rules
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## rownames: footprint length
  ## colnames: frame (0, 1, 2)
  offsets <- read.table(offsets_fname, header=T)
  offsets <- data.frame(frame=as.vector(mapply(rep, 0:2, nrow(offsets))),
                        length=rep(as.numeric(rownames(offsets)), 3),
                        offset=c(offsets$frame_0, offsets$frame_1, offsets$frame_2))
  return(offsets)
}

read_fasta_as_codons <- function(transcript_fa_fname, transcript_length_fname) {
  # read in .fasta file, convert to list of vectors of codons per transcript
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  transcript_seq <- load_fasta(transcript_fa_fname)
  transcript_lengths <- load_lengths(transcript_length_fname)
  codon_sequences <- sapply(seq_along(transcript_seq),
                            function(x) {
                              tmp_seq <- transcript_seq[x]
                              utr5_length <- subset(transcript_lengths,
                                                    transcript == names(transcript_seq)[x])$utr5_length
                              num_codons <- floor(nchar(tmp_seq) / 3)
                              sequence_offset <- utr5_length %% 3
                              codon_sequence <- substring(tmp_seq,
                                                          first=3*(seq(num_codons)-1)+1+sequence_offset,
                                                          last=3*seq(num_codons)+sequence_offset)
                              names(codon_sequence) <- as.character(seq.int(from=-floor(utr5_length/3),
                                                                            length.out=num_codons)) # start codon is 0
                              return(codon_sequence)
                            })
  names(codon_sequences) <- names(transcript_seq)
  return(codon_sequences)
}

get_codons <- function(transcript_name, cod_idx, utr5_length, transcript_seq) {
  # return codons corresponding to A-, P-, and E-site
  # for footprint originating from transcript_name and A site codon cod_idx
  ## transcript_name: character; correspond to names(transcript_seq)
  ## cod_idx: integer; codon index for A site codon
  ## utr5_length: integer; length of 5' UTR
  ## transcript_seq: character vector; transcript (+ 5' and 3' UTR regions) sequences
  A_start <- utr5_length + 3*(cod_idx-1) + 1
  A_end <- utr5_length + 3*cod_idx
  A_codon <- substr(transcript_seq[transcript_name], A_start, A_end)
  P_start <- utr5_length + 3*(cod_idx-2) + 1
  P_end <- utr5_length + 3*(cod_idx-1)
  P_codon <- substr(transcript_seq[transcript_name], P_start, P_end)
  E_start <- utr5_length + 3*(cod_idx-3) + 1
  E_end <- utr5_length + 3*(cod_idx-2)
  E_codon <- substr(transcript_seq[transcript_name], E_start, E_end)
  codons <- c(A_codon, P_codon, E_codon)
  names(codons) <- c("A", "P", "E")
  return(codons)
}

get_bias_seq <- function(dat, transcript_seq, bias_region, bias_length=3,
                         read_type="monosome") {
  ## dat: data.frame; contains columns c("transcript", "d5", "d3", "utr5_length")
  ### - if read_type=="monosome", then A site codon indices are "cod_idx"
  ### - if read_type=="disome", then A site codon indices are "cod_idx_lagging" and "cod_idx_leading"
  ## transcript_seq: character; output from load_fasta()
  ## bias_region: character vector; one of "f5" or "f3" (corresponding to 5' or 3' bias sequence)
  ## bias_length: integer; length of bias sequence
  ## read_type: character; one of "monosome" or "disome"
  ### TODO: check that read_type is valid
  if(bias_region=="f5") {
    if(read_type=="monosome") {
      seq_start <- with(dat, utr5_length + 3*(cod_idx-1)+1 - d5)
    } else {
      seq_start <- with(dat, utr5_length + 3*(cod_idx_lagging-1)+1 - d5)
    }
    seq_end <- seq_start + bias_length - 1
  } else {
    if(read_type=="monosome") {
      seq_end <- with(dat, utr5_length + 3*cod_idx + d3)
    } else {
      seq_end <- with(dat, utr5_length + 3*cod_idx_leading + d3)
    }
    seq_start <- seq_end - bias_length + 1
  }
  bias_seq <- mapply(substr, transcript_seq[as.character(dat$transcript)],
                     seq_start, seq_end)
  return(bias_seq)
}

convert_to_aa <- function(transcript_fa_fname, transcript_length_fname, output_fname) {
  # convert transcript sequence from DNA to amino acid
  ## transcript_fa_fname: character; file.path to transcript fasta (DNA) file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## output_fname: character; file.path to write output fasta amino acid sequences
  # 1. load in transcript DNA sequence
  transcript_seq <- load_fasta(transcript_fa_fname)
  transcript_lengths <- load_lengths(transcript_length_fname)
  dna_seq <- sapply(seq_along(transcript_seq),
                    function(x) {
                      tmp_seq <- transcript_seq[x]
                      tmp_name <- names(transcript_seq)[x]
                      tmp_utr5 <- subset(transcript_lengths, transcript==tmp_name)$utr5_length
                      tmp_utr3 <- subset(transcript_lengths, transcript==tmp_name)$utr3_length
                      cds_seq <- substr(tmp_seq, tmp_utr5+1, nchar(tmp_seq)-tmp_utr3)
                      return(cds_seq)
                    })
  names(dna_seq) <- names(transcript_seq)
  dna_seq <- Biostrings::DNAStringSet(dna_seq)
  # 2. translate to amino acid sequence
  aa_seq <- Biostrings::translate(dna_seq)
  # 3. write amino acid sequence to fasta file
  Biostrings::writeXStringSet(aa_seq, output_fname)
}

match_rows <- function(x, y, x_col=NULL, y_col=NULL) {
  # return row numbers of (first) matches of x in y
  ## x: data.frame
  ## y: data.frame
  ## x_col: character vector of columns to be used; NULL to use all columns
  ## y_col: character vector of columns corresponding to x_col; NULL to use x_col
  if(is.null(x_col)) {
    x_col <- colnames(x)
  }
  if(is.null(y_col)) {
    y_col <- x_col
  }
  x_terms <- do.call(paste, c(x[, x_col], sep="\r"))
  y_terms <- do.call(paste, c(y[, y_col], sep="\r"))
  match(x_terms, y_terms)
}

parse_coefs <- function(nb_fit) {
  # parse data.frame of regression coefficients for downstream analyses
  ## nb_fit: negbin object from running glm.nb()
  regression_terms <- attributes(nb_fit$terms)$term.labels
  # 1. pull coefficients from nb_fit object
  fit_coefs <- data.frame(rownames(summary(nb_fit)$coefficient),
                          summary(nb_fit)$coefficient,
                          row.names=NULL)
  colnames(fit_coefs) <- c("name", "estimate", "std_error", "t", "p")
  # 2. add reference levels
  fit_coefs <- rbind(fit_coefs,
                     data.frame(name=sapply(seq_along(names(nb_fit$xlevels)),
                                            function(x) {
                                              paste0(names(nb_fit$xlevels)[x], nb_fit$xlevels[[x]][1])
                                            }),
                                estimate=0, std_error=NA, t=NA, p=NA,
                                row.names=NULL))
  # 3. add annotations for individual terms
  for(tmp_coef in grep(":", regression_terms, invert=T, value=T)) {
    tmp_index <- grepl(paste0("^", tmp_coef), fit_coefs$name) & !grepl(":", fit_coefs$name)
    fit_coefs$group[tmp_index] <- tmp_coef
    fit_coefs$term[tmp_index] <- sub(paste0("^", tmp_coef), "", fit_coefs$name[tmp_index])
  }
  # 4. add annotations for interaction terms
  for(tmp_coef in grep(":", regression_terms, value=T)) {
    term_1 <- sub(":.*$", "", tmp_coef)
    term_2 <- sub("^.*:", "", tmp_coef)
    tmp_index <- grepl(paste0(term_1, ".*:", term_2), fit_coefs$name)
    fit_coefs$group[tmp_index] <- tmp_coef
    # transcript:(dataset) interaction terms
    ## term: transcript
    ## group_1: "transcript"
    ## group_2: (dataset)
    if("transcript" %in% c(term_1, term_2)) {
      if(term_1 == "transcript") { # first term is transcript
        fit_coefs$term[tmp_index] <- sub("transcript", "", sub(":.*", "", fit_coefs$name[tmp_index]))
      } else { # first term is (dataset)
        fit_coefs$term[tmp_index] <- sub(".*:transcript", "", fit_coefs$name[tmp_index])
      }
      fit_coefs$group_1[tmp_index] <- "transcript"
      fit_coefs$group_2[tmp_index] <- ifelse(term_1 == "transcript", term_2, term_1)
    }
    # A/P/E:(dataset) interaction terms
    ## term column: codon
    ## group_1: A/P/E
    ## group_2: (dataset)
    if((term_1 %in% c("A", "P", "E")) | (term_2 %in% c("A", "P", "E"))) {
      if(term_1 %in% c("A", "P", "E")) { # first term is codon site
        fit_coefs$term[tmp_index] <- substr(sub(":.*$", "", fit_coefs$name[tmp_index]), 2, 4)
        fit_coefs$group_1[tmp_index] <- term_1
        fit_coefs$group_2[tmp_index] <- term_2
      } else { # first term is (dataset)
        fit_coefs$term[tmp_index] <- substr(sub("^.*:", "", fit_coefs$name[tmp_index]), 2, 4)
        fit_coefs$group_1[tmp_index] <- term_2
        fit_coefs$group_2[tmp_index] <- term_1
      }
    }
    # f5:d5 / f3:d3 interaction terms
    ## term column: bias sequence
    ## group_1: f5/f3
    if(grepl("f(3|5)", tmp_coef) & grepl("d(3|5)", tmp_coef)) {
      if(grepl("f", term_1)) { # first term is bias sequence
        fit_coefs$term[tmp_index] <- sub("^.*f(3|5)", "", sub(":.*$", "", fit_coefs$name[tmp_index]))
        fit_coefs$group_1[tmp_index] <- term_1
        fit_coefs$group_2[tmp_index] <- sub(term_2, paste0(term_2, "="),
                                            sub("^.*:", "", fit_coefs$name[tmp_index]))
      } else { # first term is digest length
        fit_coefs$term[tmp_index] <- sub(".*:.*f(3|5)", "", fit_coefs$name[tmp_index])
        fit_coefs$group_1[tmp_index] <- term_2
        fit_coefs$group_2[tmp_index] <- sub(term_1, paste0(term_1, "="),
                                            sub(":.*$", "", fit_coefs$name[tmp_index]))
      }
    }
  }
  fit_coefs$group <- sub("genome_", "", fit_coefs$group)
  fit_coefs$group_1 <- sub("genome_", "", fit_coefs$group_1)
  return(fit_coefs)
}

compute_rpf_gc <- function(dat, omit="APE", transcript_fa_fname,
                           transcript_lengths_fname, read_type="monosome") {
  # compute GC content in RPF, omitting A/P/E sites
  # dat: data.frame; contains columns c("transcript", "cod_idx", "d5", "d3")
  ## omit: character; which codon sites to omit, one of c("A", "AP", "APE")
  ## transcript_fa_fname: character; file path to transcriptome fasta file
  ## transcript_lengths_fname: character; file path to transcript lengths file
  ## read_type: character; one of "monosome" or "disome"
  ### TODO: check whether codon index label matches read_type
  transcript_seq <- load_fasta(transcript_fa_fname)
  transcript_lengths <- load_lengths(transcript_lengths_fname)
  utr5_lengths <- transcript_lengths$utr5_length
  names(utr5_lengths) <- transcript_lengths$transcript
  num_omit_codons <- ifelse(grepl("E", omit), 2,
                            ifelse(grepl("P", omit), 1, 0))
  dat$transcript <- as.character(dat$transcript)
  dat$d5 <- as.numeric(as.character(dat$d5))
  dat$d3 <- as.numeric(as.character(dat$d3))
  if(read_type=="monosome") {
    A_start <- utr5_lengths[dat$transcript] + 1 + 3*(dat$cod_idx-1)
    rpf_length <- with(dat, d5+d3+3)
    # 1. extract RPF 5' regions
    rpf_5_start <- A_start - dat$d5
    rpf_5_end <- A_start - 1 - 3*num_omit_codons
    rpf_5 <- mapply(substr, transcript_seq[dat$transcript],
                    rpf_5_start, rpf_5_end)
    # 2. extract RPF 3' regions
    rpf_3_start <- A_start + 3
    rpf_3_end <- A_start + 2 + dat$d3
    rpf_3 <- mapply(substr, transcript_seq[dat$transcript],
                    rpf_3_start, rpf_3_end)
    # 3. concatenate regions
    rpf_regions <- paste0(rpf_5, rpf_3)
  } else {
    A_lagging_start <- utr5_lengths[dat$transcript] + 1 +
      3*(dat$cod_idx_lagging-1)
    A_leading_start <- utr5_lengths[dat$transcript] + 1 +
      3*(dat$cod_idx_leading-1)
    rpf_length <- with(dat, d5+d3+3*(cod_idx_leading-cod_idx_lagging+1))
    # 1. extract region 5' of lagging A site
    rpf_5_start <- A_lagging_start - dat$d5
    rpf_5_end <- A_lagging_start - 1 - 3*num_omit_codons
    rpf_5 <- mapply(substr, transcript_seq[dat$transcript],
                    rpf_5_start, rpf_5_end)
    # 2. extract region between lagging and leading A sites
    rpf_inner_start <- A_lagging_start + 3
    rpf_inner_end <- A_leading_start - 1 - 3*num_omit_codons
    rpf_inner <- mapply(substr, transcript_seq[dat$transcript],
                               rpf_inner_start, rpf_inner_end)
    # 3. extract RPF 3' regions
    rpf_3_start <- A_leading_start + 3
    rpf_3_end <- A_leading_start + 2 + dat$d3
    rpf_3 <- mapply(substr, transcript_seq[dat$transcript],
                    rpf_3_start, rpf_3_end)
    # 4. concatenate regions
    rpf_regions <- paste0(rpf_5, rpf_inner, rpf_3)
  }
  # compute GC content
  gc_content <- strsplit(rpf_regions, split="")
  gc_content <- sapply(gc_content, function(x) sum(x %in% c("G", "C")))
  gc_content <- gc_content / rpf_length
  return(gc_content)
}
