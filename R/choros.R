choros <- function(bam_fname, transcript_fa_fname, transcript_length_fname,
                   offsets_fname, num_genes=150,
                   min_prop=0.9, f5_length=2, f3_length=3,
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
  transcript_lengths <- load_lengths(transcript_length_fname)
  # 1. read in footprint alignments
  bam_data <- load_bam(bam_fname, transcript_fa_fname, transcript_length_fname,
                      offsets_fname, f5_length, f3_length)
  # 2. compute size/frame subsets
  d5_d3 <- count_d3_f3(bam_data)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3",
                                 d5_d3_subsets[x], sep="_")
                         })
  # 3. choose training set: highest RPF density (RPF count / # aa)
  ## TODO: write function to do this
  transcript_counts <- aggregate(count ~ transcript, data=bam_data, FUN=sum)
  transcript_counts$length_aa <- transcript_lengths$cds_length[match(transcript_counts$transcript,
                                                                     transcript_lengths$transcript)] / 3
  transcript_counts$TE <- with(transcript_counts, count / length_aa)
  transcript_counts <- transcript_counts[order(transcript_counts$TE, decreasing=T),]
  training_set <- as.character(transcript_counts$transcript)[1:num_genes]
  # 4. initialize data frame for regression
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets, f5_length=f5_length,
                             f3_length=f3_length, which_transcripts=training_set)
  training_data$transcript <- relevel(tunney_training$transcript, ref=training_set[1])
  training_data$count <- count_footprints(bam_data, training_data, "count")
  # 5. compute regression
  nb_fit <- MASS::glm.nb(model, data=training_data)
  # 6. correct counts
  bam_data$corrected_count <- correct_bias(bam_data, nb_fit)
  training_data$corrected_count <- count_footprints(bam_data, training_data, "corrected_count")
  # 7. evaluate and plot bias
  codon_corr <- evaluate_bias(training_data, which_column="corrected_count",
                              transcript_fa_fname, transcript_length_fname,
                              type="codon")
  bias_plot <- plot_bias(codon_corr, type="codon")
  # 8. return outputs
  output <- list(bam_alignment = bam_data,
                 model_fit = nb_fit,
                 codon_bias = codon_corr,
                 bias_plot = bias_plot)
  return(output)
}
