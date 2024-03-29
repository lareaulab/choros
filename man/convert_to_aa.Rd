% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{convert_to_aa}
\alias{convert_to_aa}
\title{Convert from DNA to amino acid sequence}
\usage{
convert_to_aa(transcript_fa_fname, transcript_length_fname, output_fname)
}
\arguments{
\item{transcript_fa_fname}{character; file path to transcriptome .fasta (DNA) file}

\item{transcript_length_fname}{character; file path to transcriptome lengths file}

\item{output_fname}{character; file path for output amino acid sequences}
}
\description{
This function extracts CDS sequence from a FASTA-formatted file, 
translates DNA into amino acid sequence, and writes amino acid sequences to disk.
}
