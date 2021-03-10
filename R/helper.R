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
  transcript_sequences <- Biostrings::readDNAStringSet(transcripts_fa_fname)
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
