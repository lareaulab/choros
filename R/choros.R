choros <- function(bam_fname, transcript_fa_fname,
                   transcript_length_fname, offsets_fname,
                   expt_dir=".", save_intermediate=T, extra_args=NULL) {
  # wrapper function to run choros pipeline
  ## bam_fname: character; file.path to .bam alignment file
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## expt_dir: character; file.path to directory to save output/intermediate files
  ## save_intermediate: logical; whether to save intermediate files
  ## extra_args : named list; additional parameters to pass onto internal functions
  ### if False: only save bam file with corrected weights
  ## TODO: test for appropriate inputs
  ## TODO: flesh out function default settings
  ## 0. parse extra_args
  default_args <- list(num_genes=250, min_prop=0.9, f5_length=3, f3_length=3,
                       model=formula(count ~ transcript + A + P + E + d5*f5 + d3*f3 + gc),
                       compute_gc=T, gc_omit="APE", num_cores=NULL, load_bam_full=F,
                       d5_d3_plot_title="", diagnostic_plot_title="")
  if(!is.null(extra_args)) {
    extra_args <- default_args
  } else {
    not_specified <- names(default_args)[!(names(default_args) %in% names(extr_args))]
    for(x in not_specified) {
      extra_args[[x]] <- default_args[[x]]
    }
  }
  # 1. generate diagnostic plot
  print("Creating diagnostic plot")
  transcript_lengths <- load_lengths(transcript_length_fname)
  diagnostic_plot <- plot_diagnostic(bam_fname, transcript_length_fname,
                                     plot_title=default_args$diagnostic_plot_title)
  # 1. read in footprint alignments
  print("Loading alignment file")
  alignment_data <- load_bam(bam_fname, transcript_fa_fname,
                             transcript_length_fname, offsets_fname,
                             f5_length=default_args$f5_length,
                             f3_length=default_args$f3_length,
                             compute_gc=default_args$compute_gc,
                             gc_omit=default_args$gc_omit,
                             num_cores=default_args$num_cores,
                             full=default_args$load_bam_full)
  # 2. compute size/frame subsets
  print("Computing d5/d3 subsets")
  d5_d3 <- count_d5_d3(alignment_data, plot_title=default_args$d5_d3_plot_title)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3",
                                 d5_d3_subsets$d3[x], sep="_")
                         })
  # 3. choose training set: highest RPF density (RPF count / # aa)
  # TODO: filter for d5/d3 subsets before calculating
  # TODO: filter for transcripts that are too short (cds_length <= exclude_codons5 + exclude_codons3 + ???)
  print("Choosing transcripts for training")
  transcript_counts <- calculate_transcript_density(alignment_data,
                                                    transcript_length_fname,
                                                    statistic=mean)
  training_set <- names(transcript_counts[1:num_genes])
  # 4. initialize data frame for regression
  print("Initializing regression data")
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=d5_d3_subsets, f5_length=f5_length,
                             f3_length=f3_length, which_transcripts=training_set)
  training_data$transcript <- relevel(training_data$transcript, ref=training_set[1])
  training_data$count <- count_footprints(alignment_data, training_data, "count")
  if(save_intermediate) {
    save(training_data, file=file.path(expt_dir,
                                       paste0(expt_name, "_training.Rda")))
  }
  # 5. compute regression
  # TODO: output warning if outlier regression values (outside -5:5 ?)
  print("Computing regression fit")
  regression_fit <- MASS::glm.nb(model, data=training_data, model=F)
  regression_coefs <- parse_coefs(regression_fit)
  # 6. correct counts
  print("Correcting RPF counts")
  alignment_data$corrected <- correct_bias(alignment_data, regression_coefs)
  training_data$corrected <- count_footprints(alignment_data, training_data,
                                              "corrected")
  # 7. evaluate bias before/after correction
  print("Evaluating bias")
  to_evaluate <- expand.grid(correction=c("count", "corrected"),
                             type=c("codon", "nt"))
  to_evaluate <- apply(to_evaluate, 1, paste, collapse="_")
  positional_bias <- sapply(to_evaluate,
                            function(x) {
                              tmp_correction <- split(x, split="_")[[1]][1]
                              tmp_type <- split(x, split="_")[[1]][2]
                              tmp_bias <- evaluate_bias(training_data,
                                                        which_column=tmp_correction,
                                                        transcript_fa_fname,
                                                        transcript_length_fname,
                                                        type=tmp_type)
                            }, simplify=F, USE.NAMES=T)
  codon_corr_max <- max(unlist(positional_bias[grepl("codon", names(positional_bias))]))
  nt_corr_max <- max(unlist(positional_bias[grepl("nt", names(positional_bias))]))
  # 8. generate plots
  codon_raw <- plot_bias(positional_bias$count_codon, type="codon") +
    coord_cartesian(ylim=c(0, codon_corr_max))
  codon_corrected <- plot_bias(positional_bias$corrected_codon, type="codon") +
    coord_cartesian(ylim=c(0, codon_corr_max))
  nt_raw <- plot_bias(positional_bias$count_nt, type="nt") +
    coord_cartesian(ylim=c(0, nt_corr_max))
  nt_corrected <- plot_bias(positional_bias$corrected_nt, type="nt") +
    coord_cartesian(ylim=c(0, nt_corr_max))
  # 9. return outputs
  all_plots <- list(diagnostic = diagnostic_plot,
                    codon_raw = codon_corr_raw_plot,
                    codon_corrected = codon_corr_corrected_plot,
                    nt_raw = nt_corr_raw_plot,
                    nt_corrected = nt_corr_corrected_plot)
  if(save_intermediate) {
    save(alignment_data,
         file=file.path(expt_dir, paste0(expt_name, "_alignments.Rda")))
    save(training_data,
         file=file.path(expt_dir, paste0(expt_name, "_training.Rda")))
    save(regression_coefs,
         file=file.path(expt_dir, paste0(expt_name, "_coefs.Rda")))
    save(positional_bias,
         file=file.path(expt_dir, paste0(expt_name, "_bias.Rda")))
    save(all_plots,
         file=file.path(expt_dir, paste0(expt_name, "_plots.Rda")))
  }
  output <- list(bam_alignment = alignment_data,
                 regression_data = training_data,
                 regression_coefs = regression_coefs,
                 plots = all_plots)
  return(output)
}
