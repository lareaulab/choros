choros <- function(bam_fname, transcript_fa_fname,
                   transcript_length_fname, offsets_fname,
                   num_genes=200, min_prop=0.9, f5_length=2, f3_length=3,
                   model=formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)) {
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
  diagnostic_plot <- plot_diagnostic(bam_fname, transcript_length_fname)
  transcript_lengths <- load_lengths(transcript_length_fname)
  # 1. read in footprint alignments
  print("Loading alignment file")
  bam_data <- load_bam(bam_fname, transcript_fa_fname, transcript_length_fname,
                       offsets_fname, f5_length, f3_length)
  # 2. compute size/frame subsets
  print("Computing d5/d3 subsets")
  d5_d3 <- count_d5_d3(bam_data)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3",
                                 d5_d3_subsets$d3[x], sep="_")
                         })
  # 3. choose training set: highest RPF density (RPF count / # aa)
  print("Choosing transcripts for training")
  transcript_counts <- calculate_transcript_density(bam_dat, transcript_length_fname,
                                                    statistic=mean)
  training_set <- names(transcript_counts[1:num_genes])
  # 4. initialize data frame for regression
  print("Initializing regression data")
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets, f5_length=f5_length,
                             f3_length=f3_length, which_transcripts=training_set)
  training_data$transcript <- relevel(tunney_training$transcript, ref=training_set[1])
  training_data$count <- count_footprints(bam_data, training_data, "count")
  # 5. compute regression
  print("Computing regression fit")
  nb_fit <- MASS::glm.nb(model, data=training_data, model=F)
  fit_coefs <- parse_coefs(nb_fit)
  # 6. correct counts
  print("Correcting RPF counts")
  bam_data$corrected_count <- correct_bias(bam_data, nb_fit)
  training_data$corrected_count <- count_footprints(bam_data, training_data, "corrected_count")
  # 7. evaluate and plot bias
  print("Evaluating bias")
  codon_corr <- lapply("count", "corrected_count",
                       function(x) {
                         evaluate_bias(training_data, which_column=x,
                                       transcript_fa_fname, transcript_length_fname,
                                       type="codon")
                       })
  codon_corr_plots <- lapply(codon_corr, plot_bias)
  nt_corr <- lapply("count", "corrected_count",
                    function(x) {
                      evaluate_bias(training_data, which_column=x,
                                    transcript_fa_fname, transcript_length_fname,
                                    type="nt")
                    })
  nt_corr_plots <- lapply(nt_corr, plot_bias, type="nt")
  # 8. return outputs
  output <- list(diagnostic_plot = diagnostic_plot,
                 bam_alignment = bam_data,
                 regression_data = training_data,
                 model_fit = nb_fit,
                 regression_coefs = fit_coefs,
                 codon_bias = codon_corr,
                 nt_bias = nt_corr,
                 codon_plots = codon_corr_plots,
                 nt_plots = nt_corr_plots)
  return(output)
}
