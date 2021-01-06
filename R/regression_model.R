model_frame_by_size <- function(bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                                min_prop=0.9, num_regression_genes=100, f5_length=3, f3_length=3, 
                                model=formula(count ~ transcript + A + P + E + f5 + f3), 
                                plot_title="", trunc5=20, trunc3=20, num_f5_codons=6, num_f3_codons=6) {
  # wrapper function
  ## bam_fname: character; file.path to .bam alignment file
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## min_prop: numeric ; minimum proportion of footprints to be represented in d5/d3 subsets
  ## num_genes: integer; number to genes to process in regression modeling
  ## f5_length: integer; length of 5' bias region
  ## f3_length: integer; length of 3' bias region
  ## model: formula; formula for MASS::glm.nb() to be applied to each data subset
  ## plot_title: character; title for output plots
  ## save_models: logical; whether to save list of regression models as .Rda file
  ## output_dir: character; file.path to output directory
  ## trunc5: integer; number of 5' codons to leave out in codon correlation regression
  ## trunc3: integer; number of 3' codons to leave out in codon correlation regression
  ## num_f5_codons: integer; number of codons to include 5' of A site in regression for codon correlation
  ## num_f3_codons: integer; number of codons to include 3' of A site in regression for codon correlation
  # 1. load .bam alignment
  footprints <- load_bam(bam_fname, transcript_length_fname, offsets_fname)
  # 2. compute proportion of footprints per d5/d3 subset
  print("Computing subsets")
  d5_d3 <- count_d5_d3(footprints)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  # 3. initiate data.frames for regression
  transcript_counts <- aggregate(count ~ transcript, data=footprints, FUN=sum)
  transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]
  top_transcripts <- as.character(transcript_counts$transcript)[1:num_regression_genes]
  print("Initiating regression data")
  regression_data <- init_data(transcript_fa_fname, transcript_length_fname, 
                               d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                               which_transcripts=top_transcripts)
  print("Counting footprints")
  regression_data$count <- count_footprints(footprints, regression_data)
  # 4. compute regression
  subset_data <- lapply(seq(nrow(d5_d3_subsets)),
                        function(x) {
                          subset(regression_data, d5==d5_d3_subsets$d5[x] & d3==d5_d3_subsets$d3[x])
                        })
  print("Computing regression models")
  subset_models <- lapply(subset_data, function(x) { MASS::glm.nb(model, data=x, model=F) })
  names(subset_models) <- subset_names
  # 5. evaluate bias correction
  print("Evaluating bias correction")
  subset_data <- do.call(rbind,
                         lapply(seq(nrow(d5_d3_subsets)),
                                function(x) {
                                  correct_bias(subset_data[[x]], subset_models[[x]])
                                }))
  transcript_lengths <- load_lengths(transcript_length_fname)
  leaveOneOut_codon_corr <- evaluate_bias(subset_data, transcript_fa_fname, 
                                          unique(transcript_lengths$utr5_length),
                                          unique(transcript_lengths$utr3_length),
                                          trunc5, trunc3, num_f5_codons, num_f3_codons, type="codon")
  leaveOneOut_nt_corr <- evaluate_bias(subset_data, transcript_fa_fname, 
                                       unique(transcript_lengths$utr5_length),
                                       unique(transcript_lengths$utr3_length),
                                       trunc5, trunc3, num_f5_codons, num_f3_codons, type="nt")
  leaveOneOut_codon_plot <- plot_bias(leaveOneOut_codon_corr, plot_title)
  leaveOneOut_nt_plot <- plot_bias(leaveOneOut_nt_corr, plot_title)
  # 7. correct footprints
  print("Correcting footprint counts")
  footprints <- correct_bias_bam(footprints, subset_models)
  # 8. return output
  return(list(footprints=footprints, regression_data=subset_data, 
              subsets=d5_d3$counts, models=subset_models, 
              subsets_plot=d5_d3$plot, 
              codon_corr_plot=leaveOneOut_codon_plot, nt_corr_plot=leaveOneOut_nt_plot))
}