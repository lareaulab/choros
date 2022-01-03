choros <- function(bam_fname, transcript_fa_fname,
                   transcript_length_fname, offsets_fname,
                   num_genes=200, min_prop=0.9, f5_length=2, f3_length=3,
                   model=formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3),
                   expt_dir=".") {
  # wrapper function to run choros pipeline
  ## bam_fname: character; file.path to .bam alignment file
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## min_prop: numeric; proportion of footprints to be modeled by d5/d3 lengths
  ## f5_length: integer; length of 5' bias region
  ## f3_length: integer; length of 3' bias region
  ## num_genes: integer; number of genes used in regression model
  ## TODO: test for appropriate inputs
  ## TODO: flesh out function default settings
  print("Creating diagnostic plot")
  transcript_lengths <- load_lengths(transcript_length_fname)
  diagnostic_plot <- plot_diagnostic(bam_fname, transcript_length_fname)
  save(diagnostic_plot, file=file.path(expt_dir, "diagnostic_plot.Rda"))
  # 1. read in footprint alignments
  print("Loading alignment file")
  alignment_data <- load_bam(bam_fname, transcript_fa_fname, transcript_length_fname,
                             offsets_fname, f5_length, f3_length)
  save(alignment_data, file=file.path(expt_dir, "alignment.Rda"))
  # 2. compute size/frame subsets
  print("Computing d5/d3 subsets")
  d5_d3 <- count_d5_d3(alignment_data)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3",
                                 d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file=file.path(expt_dir, "d5_d3_subsets.Rda"))
  # 3. choose training set: highest RPF density (RPF count / # aa)
  # TODO: filter for d5/d3 subsets before calculating
  # TODO: filter for transcripts that are too short (cds_length <= exclude_codons5 + exclude_codons3 + ???)
  print("Choosing transcripts for training")
  transcript_counts <- calculate_transcript_density(alignment_data, transcript_length_fname,
                                                    statistic=mean)
  training_set <- names(transcript_counts[1:num_genes])
  writeLines(training_set, file.path(expt_dir, "training_transcripts.txt"))
  save(transcript_counts, file=file.path(expt_dir, "transcript_counts"))
  # 4. initialize data frame for regression
  print("Initializing regression data")
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=d5_d3_subsets, f5_length=f5_length,
                             f3_length=f3_length, which_transcripts=training_set)
  training_data$transcript <- relevel(training_data$transcript, ref=training_set[1])
  training_data$count <- count_footprints(alignment_data, training_data, "count")
  save(training_data, file=file.path(expt_dir, "training_data.Rda"))
  # 5. evaluate bias before correction
  codon_corr_raw <- evaluate_bias(training_data, "count",
                                  transcript_fa_fname, transcript_length_fname,
                                  type="codon")
  codon_corr_raw_plot <- plot_bias(codon_corr_raw)
  save(codon_corr_raw, codon_corr_raw_plot,
       file=file.path(expt_dir, "codon_corr_raw.Rda"))
  nt_corr_raw <- evaluate_bias(training_data, "count",
                               transcript_fa_fname, transcript_length_fname,
                               type="nt")
  nt_corr_raw_plot <- plot_bias(nt_corr_raw)
  save(nt_corr_raw, nt_corr_raw_plot,
       file=file.path(expt_dir, "nt_corr_raw.Rda"))
  # 5. compute regression
  print("Computing regression fit")
  regression_fit <- MASS::glm.nb(model, data=training_data, model=F)
  regression_coefs <- parse_coefs(regression_fit)
  # TODO: output warning if outlier regression values (outside -5:5 ?)
  save(regression_fit, file=file.path(expt_dir, "regression_fit.Rda"))
  save(regression_coefs, file=file.path(expt_dir, "regression_coefs.Rda"))
  # 6. correct counts
  print("Correcting RPF counts")
  alignment_data$corrected_count <- correct_bias(alignment_data, regression_fit)
  training_data$corrected_count <- count_footprints(alignment_data, training_data, "corrected_count")
  save(alignment_data, file=file.path(expt_dir, "alignment_data.Rda"))
  save(training_data, file=file.path(expt_dir, "training_data.Rda"))
  # 5. evaluate bias before correction
  codon_corr_corrected <- evaluate_bias(training_data, "corrected_count",
                                        transcript_fa_fname, transcript_length_fname,
                                        type="codon")
  codon_corr_corrected_plot <- plot_bias(codon_corr_corrected)
  save(codon_corr_corrected, codon_corr_corrected_plot,
       file=file.path(expt_dir, "codon_corr_corrected.Rda"))
  nt_corr_corrected <- evaluate_bias(training_data, "corrected_count",
                                     transcript_fa_fname, transcript_length_fname,
                                     type="nt")
  nt_corr_corrected_plot <- plot_bias(nt_corr_corrected)
  save(nt_corr_corrected, nt_corr_corrected_plot,
       file=file.path(expt_dir, "nt_corr_corrected.Rda"))
  # 8. return outputs
  output <- list(diagnostic_plot = diagnostic_plot,
                 bam_alignment = alignment_data,
                 regression_data = training_data,
                 model_fit = regression_fit,
                 regression_coefs = regression_coefs,
                 bias = list(codon_raw = codon_corr_raw,
                             codon_corrected = codon_corr_corrected,
                             nt_raw = nt_corr_raw,
                             nt_corrected = nt_corr_corrected),
                 plots = list(codon_raw = codon_corr_raw_plot,
                              codon_corrected = codon_corr_corrected_plot,
                              nt_raw = nt_corr_raw_plot,
                              nt_corrected = nt_corr_corrected_plot))
  return(output)
}
