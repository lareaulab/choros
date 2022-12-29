#' Generate diagnostic plots
#' 
#' @description 
#' This function generates diagnostic plots: RPF density at the start
#' codon; RPF density at the stop codon; histogram of RPF lengths; heatmap of 
#' RPF counts by RPF length and frame
#' 
#' Input `.bam` file must contain `tag.ZW` of posterior weights computed from
#' `RSEM`
#' 
#' @param bam_fname character; file path to .bam alignment file
#' @param transcript_length_fname character; file path to transcriptome lengths file
#' @param plot_title character; title for diagnostic plot
#' @param start_min integer; 5'-most position to include in plot, relative to start codon
#' @param start_max integer; 3'-most position to include in plot, relative to start codon
#' @param stop_min integer; 5'-most position to include in plot, relative to stop codon
#' @param stop_max integer; 3'-most position to include in plot, relative to stop codon
#' @param length_min integer; minimum RPF length
#' @param length_max integer; maximum RPF length
#' 
#' @returns A patchwork object of diagnostic plots
plot_diagnostic <- function(bam_fname, transcript_length_fname, plot_title="",
                            start_min=-20, start_max=5,
                            stop_min=-5, stop_max=20,
                            length_min=15, length_max=35) {
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
  bam_dat$frame <- factor(with(bam_dat, (pos - utr5_length - 1) %% 3),
                          levels=0:2)
  bam_dat$start_distance <- with(bam_dat, pos-(utr5_length+1+3))
  bam_dat$stop_distance <- with(bam_dat, (pos+qwidth-1)-(utr5_length+cds_length))
  num_reads <- round(sum(bam_dat$tag.ZW))
  # start codon plot
  start_ribogrid <- subset(bam_dat, (start_distance %in% seq(start_min, start_max)) &
                             (qwidth %in% seq(length_min, length_max)))
  start_ribogrid <- aggregate(tag.ZW ~ start_distance + qwidth,
                              data=start_ribogrid, FUN=sum)
  plot_start_ribogrid <- ggplot(start_ribogrid,
                                aes(x=start_distance, y=qwidth, fill=tag.ZW)) +
    geom_tile(color=1) + theme_classic() +
    scale_fill_gradient(low="white", high="blue", name="count") +
    xlab("read 5' end to start codon (P site)") + ylab("RPF length") +
    scale_x_continuous(breaks=seq(from=start_min, to=start_max, by=5)) +
    scale_y_continuous(breaks=seq(from=length_min, to=length_max, by=5))
  start_hist_dat <- aggregate(tag.ZW ~ start_distance, data=start_ribogrid, FUN=sum)
  start_hist_dat$tag.ZW <- with(start_hist_dat, tag.ZW / sum(tag.ZW))
  plot_start <- ggplot(start_hist_dat, aes(x=start_distance, y=tag.ZW)) +
    geom_col() + theme_classic() +
    xlab("distance between read 5' end and start codon") + ylab("density")
  # stop codon plot
  stop_ribogrid <- subset(bam_dat, (stop_distance %in% seq(stop_min, stop_max)) &
                            (qwidth %in% seq(length_min, length_max)))
  stop_ribogrid <- aggregate(tag.ZW ~ stop_distance + qwidth,
                             data=stop_ribogrid, FUN=sum)
  plot_stop_ribogrid <- ggplot(stop_ribogrid,
                               aes(x=stop_distance, y=qwidth, fill=tag.ZW)) +
    geom_tile(color=1) + theme_classic() +
    scale_fill_gradient(low="white", high="blue", name="count") +
    xlab("read 3' end to stop codon") + ylab("RPF length") +
    scale_x_continuous(breaks=seq(from=stop_min, to=stop_max, by=5)) +
    scale_y_continuous(breaks=seq(from=length_min, to=length_max, by=5))
  stop_hist_dat <- aggregate(tag.ZW ~ stop_distance, data=stop_ribogrid, FUN=sum)
  stop_hist_dat$tag.ZW <- with(stop_hist_dat, tag.ZW / sum(tag.ZW))
  plot_stop <- ggplot(stop_hist_dat, aes(x=stop_distance, y=tag.ZW)) +
    geom_col() + theme_classic() +
    xlab("distance between read 3' end and stop codon") + ylab("density")
  # histogram of fragment lengths
  length_dat <- aggregate(tag.ZW ~ qwidth + frame, data=bam_dat, FUN=sum)
  length_dat <- subset(length_dat, qwidth >= length_min & qwidth <= length_max)
  plot_length <- ggplot(length_dat, aes(x=qwidth, y=tag.ZW, fill=frame)) +
    geom_col() + theme_classic() + xlab("fragment length (nt)") + ylab("") +
    ggtitle(plot_title, subtitle=paste("n =", num_reads, "footprints")) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Set1")) +
    theme(legend.position="right")
  # frame v. length plot
  frame_length_dat <- aggregate(tag.ZW ~ qwidth + frame, data=bam_dat, FUN=sum)
  plot_frame_length <- ggplot(frame_length_dat, aes(x=frame, y=qwidth, fill=tag.ZW)) +
    geom_tile(color=1) + theme_classic() + coord_fixed(ratio=1) +
    scale_fill_gradient(low="white", high="blue", name="count") +
    xlab("frame") + ylab("RPF length")
  # return plots
  aggregate_plot <- plot_length +
    plot_start + plot_stop +
    plot_start_ribogrid + plot_stop_ribogrid +
    plot_frame_length +
    plot_layout(design="
                11116
                22336
                44556")
  return(aggregate_plot)
}

#' Plot RPF biases
#' 
#' @description 
#' This function generates `iXnos`-style plot of codon or nucleotide positional 
#' contributions to RPF counts.
#' 
#' @param model_metric numeric vector; output from `evaluate_bias`
#' @param plot_title character; title for plot
#' @param plot_subtitle character; subtitle for plot
#' @param type character; one of `codon` or `nt` corresponding to `model_metric`
#' @param metric character; one of `corr` or `norm`
#' @param fill_colors named character vector; colors corresponding to positions
#' 
#' @returns A ggplot object
plot_bias <- function(model_metric, plot_title="", plot_subtitle="",
                      type="codon", metric="corr", fill_colors=NULL) {
  if(is.null(fill_colors)) {
    fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")
    names(fill_colors) <- c("bias", "E", "P", "A", "other")
  } else {
    # TODO: add checks for whether fill_colors is valid
    if(length(fill_colors) == 1 & fill_colors=="none") {
      fill_colors <- rep("grey", 5)
      names(fill_colors) <- c("bias", "E", "P", "A", "other")
    }
  }
  names(model_metric) <- sub("n", "-", names(model_metric))
  names(model_metric) <- sub("p", "", names(model_metric))
  model_metric <- data.frame(position=factor(names(model_metric),
                                             levels=names(model_metric)),
                             value=model_metric,
                             col="other", stringsAsFactors=F)
  model_metric$col[grepl("^A", as.character(model_metric$position))] <- "A"
  model_metric$col[grepl("^P", as.character(model_metric$position))] <- "P"
  model_metric$col[grepl("^E", as.character(model_metric$position))] <- "E"
  if(type=="codon") {
    model_metric$col[model_metric$position %in% c(-4, -5, 3, 4)] <- "bias"
  } else {
    suppressWarnings(tmp_position <- as.numeric(as.character(model_metric$position)))
    model_metric$col[tmp_position <= -14 | tmp_position >= 9] <- "bias"
  }
  bias_plot <- ggplot(model_metric, aes(x=position, y=value, fill=col)) +
    geom_col() +theme_bw() + xlab("position") +
    ggtitle(plot_title, subtitle=plot_subtitle) + theme(legend.position="none") +
    scale_fill_manual(values=fill_colors)
  if(type=="nt") {
    bias_plot <- bias_plot + theme(axis.text.x=element_text(angle=90, vjust=0.5))
  }
  if(metric=="corr") {
    bias_plot <- bias_plot + ylab(expression(paste(Delta, " correlation")))
  } else {
    bias_plot <- bias_plot + ylab(expression(paste(Sigma, "(", beta^2, ")")))
  }
  return(bias_plot)
}

#' Visualize regression coefficients
#' 
#' @description 
#' This function generates either a violin plot of regression coefficient values
#' or a volcano plot of regression coefficient values by p-value.
#' 
#' @param model_coefs data frame; output from `parse_coefs`
#' @param plot_type character; one of `violin` or `volcano`
#' @param volcano_term character; coefficient group to plot
#' @param conf_int numeric; threshold for statistical significance
#' @param p_adj_method character; `method` argument for `p.adjust` for multiple hypothesis correction
#' @param p_adj_group character; one of `all` or `subset_only`
plot_coefs <- function(model_coefs, plot_type="violin",
                       volcano_term="A", conf_int=0.95,
                       p_adj_method="fdr", p_adj_group=c("all", "subset_only")) {
  if(plot_type=="violin") {
    model_coefs <- subset(model_coefs, !is.na(model_coefs$group))
    model_coefs$group <- factor(model_coefs$group,
                                levels=unique(model_coefs$group))
    return(ggplot(subset(model_coefs, !is.na(model_coefs$group)),
                  aes(x=group, y=estimate, fill=group)) +
             geom_hline(yintercept=0) + geom_violin() +
             theme_classic() + guides(fill="none") +
             xlab("") + ylab(expr(beta)) +
             theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))
  }
  if(plot_type=="volcano") {
    if(p_adj_group == "subset_only") {
      model_coefs <- subset(model_coefs, group_1==volcano_term)
    }
    model_coefs$p_adj <- p.adjust(model_coefs$p, method=p_adj_method)
    model_coefs$p_signif <- ifelse(model_coefs$p_adj<(1-conf_int),
                                   paste("p_adj <", 1-conf_int),
                                   paste("p_adj >=", 1-conf_int))
    color_legend <- c("red", "black")
    names(color_legend) <- c(paste("p_adj <", 1-conf_int),
                             paste("p_adj >=", 1-conf_int))
    return(ggplot(subset(model_coefs, group_1==volcano_term),
                  aes(x=estimate, y=-log10(p),
                      xmin=estimate+qnorm((1-conf_int)/2)*std_error,
                      xmax=estimate+qnorm(1-(1-conf_int)/2)*std_error,
                      col=p_signif)) +
             geom_point() + geom_errorbarh(alpha=0.5) +
             geom_hline(yintercept=0) + geom_vline(xintercept=0) +
             theme_classic() + labs(col="") + scale_color_manual(values=color_legend) +
             xlab(expr(beta)) + ylab("-log10(p)"))
  }
}
