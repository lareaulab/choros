## requires libraries: Biostrings

load_lengths <- function(lengths_fname) {
  # load transcript lengths table
  ## length_fname: character; file path to transcript lengths file
  transcript_lengths <- read.table(lengths_fname, stringsAsFactors=F,
                                   col.names=c("transcript", "utr5_length", "cds_length", "utr3_length"))
  return(transcript_lengths)
}

load_fa <- function(transcript_fa_fname) {
  # load transcript sequences from genome .fa file
  ## transcripts_fa_fname: character; file path to transcriptome .fa file
  raw_text <- readLines(transcript_fa_fname)
  transcript_startLines <- grep(">", raw_text)
  num_transcripts <- length(transcript_startLines)
  transcript_names <- sapply(transcript_startLines,
                             function(x) {
                               gsub(">", "", strsplit(raw_text[x], split=" ")[[1]][1])
                             })
  transcript_startLines <- c(transcript_startLines, length(raw_text)+1) # add extra line for bookkeeping
  transcript_sequences <- sapply(1:num_transcripts,
                                 function(x) {
                                   startLine <- transcript_startLines[x]+1
                                   endLine <- transcript_startLines[x+1]-1
                                   transcriptSequence <- paste(raw_text[startLine:endLine], collapse="")
                                   return(transcriptSequence)
                                 })
  names(transcript_sequences) <- transcript_names
  return(transcript_sequences)
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

readFAfile <- function(faFile, utr5_length, utr3_length) {
  ## read in .fa file where sequence broken up over multiple lines, convert to list of vectors of codons per transcript
  # faFile: character; path to .fa file of transcript sequences
  # utr5_length: numeric; length of 5' UTR
  # utr3_length: numeric; length of 3' UTR
  rawFile <- readLines(faFile)
  transcriptStartLines <- grep(">", rawFile)
  nTranscripts <- length(transcriptStartLines)
  transcriptNames <- sapply(transcriptStartLines,
                            function(x) {
                              gsub(">", "", strsplit(rawFile[x], split=" ")[[1]][1])
                            })
  transcriptStartLines <- c(transcriptStartLines, length(rawFile)+1, length(rawFile)+1) # add extra line for bookkeeping
  faList <- sapply(1:nTranscripts,
                   function(x) {
                     startLine <- transcriptStartLines[x]+1
                     endLine <- transcriptStartLines[x+1]-1
                     transcriptSequence <- paste(rawFile[startLine:endLine], collapse="")
                     nCodons <- floor((nchar(transcriptSequence))/3)
                     sequenceOffset <- utr5_length %% 3
                     codonSequence <- substring(transcriptSequence,
                                                first=(3*(1:nCodons-1)+1)+sequenceOffset,
                                                last=(3*(1:nCodons))+sequenceOffset)
                     names(codonSequence) <- as.character(seq.int(from=-floor(utr5_length/3), length.out=nCodons)) # start codon is 0
                     return(codonSequence)
                   })
  names(faList) <- transcriptNames
  return(faList)
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

convert_to_aa <- function(transcript_fa_fname, output_fname, utr5_length=20, utr3_length=20) {
  # convert transcript sequence from DNA to amino acid
  ## transcript_fa_fname: character; file.path to transcript fasta (DNA) file
  ## output_fname: character; file.path to write output fasta amino acid sequences
  ## utr5_length: integer; length of 5' UTR
  ## utr3_length: integer; length of 3' UTR
  # 1. load in transcript DNA sequence
  dna_seq <- load_fa(transcript_fa_fname)
  dna_seq <- sapply(dna_seq, function(x) { substr(x, utr5_length+1, nchar(x)-utr3_length) })
  dna_seq <- Biostrings::DNAStringSet(dna_seq)
  # 2. translate to amino acid sequence
  aa_seq <- Biostrings::translate(dna_seq)
  # 3. write amino acid sequence to fasta file
  output_aa <- rep(NA, 2*length(aa_seq))
  output_aa[2*(1:length(aa_seq))-1] <- paste0(">", names(aa_seq))
  output_aa[2*(1:length(aa_seq))] <- sub("*", "", as.character( aa_seq))
  writeLines(output_aa, con=output_fname)
}
