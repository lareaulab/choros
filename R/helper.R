load_lengths <- function(lengths_fname) {
  # load transcript lengths table
  ## length_fname: character; file path to transcript lengths file
  transcript_lengths <- read.table(lengths_fname, stringsAsFactors=F,
                                   col.names=c("transcript", "utr5_length", "cds_length", "utr3_length"))
  return(transcript_lengths)
}

load_gff <- function(gff_fname) {
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
                                   row.names=NULL)
  return(transcript_lengths)
}

load_fa <- function(transcript_fa_fname) {
  # load transcript sequences from genome .fa file
  ## transcripts_fa_fname: character; file path to transcriptome .fa file
  transcript_sequences <- Biostrings::readDNAStringSet(transcript_fa_fname)
  transcript_sequences <- as.character(transcript_sequences)
  return(transcript_sequences)
}

write_fa <- function(transcript_seq, fa_fname) {
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

read_fa_codons <- function(transcript_fa_fname, transcript_length_fname) {
  # read in .fa file, convert to list of vectors of codons per transcript
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  transcript_seq <- load_fa(transcript_fa_fname)
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

get_bias_seq <- function(transcript_name, cod_idx, digest_length, utr5_length,
                         transcript_seq, bias_region, bias_length=2) {
  # get bias sequence at end of footprint
  ## transcript_name: character; transcript name, corresponds with item in names(transcript_seq)
  ## cod_idx: integer; index of A site codon
  ## digest_length: integer; d5 or d3 length between A site and footprint end
  ## utr5_length: integer; length of 5' UTR region (from lengths file)
  ## transcript_seq: character vector; transcript sequences (+ 5' and 3' UTR regions)
  ## bias_region: character; f5 or f3 (corresponding to 5' or 3' bias sequence)
  ## bias_length: integer; length of bias sequence
  if(bias_region=="f5") {
    seq_start <- utr5_length + 3*(cod_idx-1)+1 - digest_length
    seq_end <- seq_start + bias_length - 1
  } else {
    if(bias_region=="f3") {
      seq_end <- utr5_length + 3*cod_idx + digest_length
      seq_start <- seq_end - bias_length + 1
    }
  }
  bias_seq <- substr(transcript_seq[transcript_name], seq_start, seq_end)
  return(bias_seq)
}

convert_to_aa <- function(transcript_fa_fname, transcript_length_fname, output_fname) {
  # convert transcript sequence from DNA to amino acid
  ## transcript_fa_fname: character; file.path to transcript fasta (DNA) file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## output_fname: character; file.path to write output fasta amino acid sequences
  # 1. load in transcript DNA sequence
  transcript_seq <- load_fa(transcript_fa_fname)
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

plot_diagnostic <- function(bam_fname, transcript_length_fname, plot_title="",
                            start_min=-20, start_max=5,
                            stop_min=-5, stop_max=20,
                            length_min=15, length_max=35) {
  # return plot of footprint density around start and stop codons
  ## bam_fname: character; file.path to bam alignment file
  ## - must contain tag.ZW from RSEM
  ## transcript_length_fname: character; file.path to transcript lengths file
  ## plot_title: character; plot title
  ## start_min: integer; 5'-most position to plot, relative to start codon
  ## start_max: integer; 3'-most position to plot, relative to start codon
  ## stop_min: integer; 5'-most position to plot, relative to stop codon
  ## stop_max: integer; 3'-most position to plot, relative to stop codon
  ## length_min: integer; minimum fragment length
  ## length_max: integer; maximum fragment length
  # load transcript lengths file
  transcript_lengths <- load_lengths(transcript_length_fname)
  # load alignment
  features <- c("rname", "pos", "seq", "qwidth")
  bam_param <- Rsamtools::ScanBamParam(tag=c("ZW", "MD"), what=features)
  bam_file <- Rsamtools::BamFile(bam_fname)
  bam_dat <- data.frame(Rsamtools::scanBam(bam_file, param=bam_param)[[1]])
  # subset to alignments to transcripts
  bam_dat <- subset(bam_dat, !is.na(rname))
  # wrangle data
  bam_dat$utr5_length <- transcript_lengths$utr5_length[match(bam_dat$rname,
                                                              transcript_lengths$transcript)]
  bam_dat$cds_length <- transcript_lengths$cds_length[match(bam_dat$rname,
                                                            transcript_lengths$transcript)]
  bam_dat$d5 <- with(bam_dat, pos - utr5_length)
  bam_dat$d3 <- with(bam_dat, (pos + qwidth - 1) - (utr5_length + cds_length))
  bam_dat$dist_stop <- with(bam_dat, pos - (utr5_length + cds_length))
  num_reads <- round(sum(bam_dat$tag.ZW))
  # start codon plot
  d5_dat <- subset(bam_dat, (d5 %in% seq(start_min, start_max)) &
                     (qwidth %in% seq(length_min, length_max)))
  d5_dat <- aggregate(tag.ZW ~ d5 + qwidth, data=d5_dat, FUN=sum)
  plot_d5 <- ggplot(d5_dat, aes(x=d5, y=qwidth, fill=tag.ZW)) +
    geom_tile(color=1) + theme_classic() +
    scale_fill_gradient(low="white", high="blue", name="count") +
    xlab("distance between read 5' end and start codon") + ylab("fragment length")
  start_dat <- subset(bam_dat, d5 > start_min & d5 < 50)
  start_dat <- aggregate(tag.ZW ~ d5, data=start_dat, FUN=sum)
  start_dat$tag.ZW <- start_dat$tag.ZW / sum(bam_dat$tag.ZW)
  plot_start <- ggplot(start_dat, aes(x=d5, y=tag.ZW)) +
    geom_col() + theme_classic() +
    xlab("distance between read 5' end and start codon") + ylab("density")
  # stop codon plot
  d3_dat <- subset(bam_dat, (d3 %in% seq(stop_min, stop_max)) &
                     (qwidth %in% seq(length_min, length_max)))
  d3_dat <- aggregate(tag.ZW ~ d3 + qwidth, data=d3_dat, FUN=sum)
  plot_d3 <- ggplot(d3_dat, aes(x=d3, y=qwidth, fill=tag.ZW)) +
    geom_tile(color=1) + theme_classic() +
    scale_fill_gradient(low="white", high="blue", name="count") +
    xlab("distance between read 3' end and stop codon") + ylab("fragment length")
  stop_dat <- subset(bam_dat, dist_stop > -60 & dist_stop < 10)
  stop_dat <- aggregate(tag.ZW ~ dist_stop, data=stop_dat, FUN=sum)
  stop_dat$tag.ZW <- stop_dat$tag.ZW / sum(bam_dat$tag.ZW)
  plot_stop <- ggplot(stop_dat, aes(x=dist_stop, y=tag.ZW)) +
    geom_col() + theme_classic() +
    xlab("distance between read 5' end and stop codon") + ylab("density")
  # histogram of fragment lengths
  length_dat <- aggregate(tag.ZW ~ qwidth, data=bam_dat, FUN=sum)
  length_dat <- subset(length_dat, qwidth >= 15 & qwidth <= 50)
  plot_length <- ggplot(length_dat, aes(x=qwidth, y=tag.ZW)) +
    geom_col() + theme_classic() + xlim(15, 50) +
    xlab("fragment length (nt)") + ylab("") +
    ggtitle(plot_title, subtitle=paste("n =", num_reads, "footprints"))
  # return plots
  aggregate_plot <- plot_length / (plot_start + plot_stop) / (plot_d5 + plot_d3)
  return(aggregate_plot)
}

calculate_codon_density <- function(bam_dat, transcript_length,
                                    exclude_codons5=10, exclude_codons3=10) {
  # return vector of counts per codon
  ## bam_dat: data.frame; output (or subset of) from load_bam()
  ## transcript_lenght: integer; number of codons in transcript
  ## exclude_codons5: integer; number of codons to exclude from 5' end of transcript
  ## exclude_codons3: integer; number of codons to exclude from 3' end of transcript
  # 1. initialize empty vector
  counts <- rep(0, transcript_length)
  names(counts) <- seq(transcript_length)
  # 2. count up footprints per codon
  bam_dat_cts <- aggregate(count ~ cod_idx, data=bam_dat, FUN=sum)
  # 3. populate counts vector
  counts[bam_dat_cts$cod_idx] <- bam_dat_cts$count
  # 4. remove first and last codons
  counts <- counts[(1+exclude_codons5):(length(counts)-exclude_codons3)]
  return(counts)
}

calculate_transcript_density <- function(bam_dat, transcript_length_fname, statistic,
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
