library('ggplot2')

args <- base::commandArgs(trailingOnly=TRUE)
comparison_tsv_path <- args[1]
out_dir <- args[2]
run_a_name <- args[3]
run_b_name <- args[4]

get_colors <- function() {
  black <- '#000000'
  blue <- '#0000FF'
  green <- '#00FF00'
  orange <- '#FFB428'
  red <- '#FF0000'
  white <- '#FFFFFF'
  return(list(black=black,
              blue=blue,
              green=green,
              orange=orange,
              red=red,
              white=white))
}

convert_column_with_none_to_numeric <- function(df, name) {
  values <- df[[name]]
  is_none <- values == 'None'
  values[is_none] <- NA
  return(base::as.numeric(values))
}

read_comparison_tsv <- function(comparison_tsv_path) {
  df <- utils::read.table(comparison_tsv_path, header=TRUE, sep='\t', quote='')
  ## id_a, id_b, sample_1_a, sample_1_b, ..., sample_n_a, sample_n_b
  headers <- base::colnames(df)
  num_headers <- base::length(headers)
  num_samples <- (num_headers - 2) / 2
  sample_names <- base::vector(mode='character')
  for (sample_i in 1:num_samples) {
    ## skip 1st two headers and then each sample has two headers
    header_i <- 3 + ((sample_i - 1) * 2)
    header_a <- headers[header_i]
    header_b <- headers[header_i + 1]
    header_len <- base::nchar(header_a)
    ## skip trailing '_a'
    sample_name <- base::substr(header_a, 1, header_len - 2)
    sample_names <- base::append(sample_names, sample_name)

    df[[header_a]] <- convert_column_with_none_to_numeric(df, header_a)
    df[[header_b]] <- convert_column_with_none_to_numeric(df, header_b)
  }

  return(base::list(df=df, sample_names=sample_names))
}

save_plot <- function(plot, name, out_dir) {
  out_path_pdf <- base::paste0(out_dir, '/', name, '_plot.pdf')
  ggplot2::ggsave(out_path_pdf, plot=plot, dpi=300, width=8, height=8)

  out_path_png <- base::paste0(out_dir, '/', name, '_plot.png')
  ggplot2::ggsave(out_path_png, plot=plot, dpi=300, width=8, height=8)
}

create_ggplot_theme <- function() {
  colors <- get_colors()
  return(
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.line=ggplot2::element_line(),
                   plot.background=ggplot2::element_rect(fill=colors$white))
  )
}

make_plot_for_sample <- function(df, sample, run_a_name, run_b_name, out_dir) {
  colors <- get_colors()
  header_a <- base::paste0(sample, '_a')
  header_b <- base::paste0(sample, '_b')
  plot_df <- base::data.frame(run_a=df[[header_a]],
                              run_b=df[[header_b]])

  ## NA abundance is converted to zero
  plot_df$run_a[base::is.na(plot_df$run_a)] <- 0
  plot_df$run_b[base::is.na(plot_df$run_b)] <- 0

  corr <- stats::cor(plot_df$run_a, plot_df$run_b, method='pearson')
  corr <- base::round(corr, digits=3)

  zoom_cutoff <- 1000
  zoom_df <- plot_df[(plot_df$run_a <= zoom_cutoff)
                     & (plot_df$run_b <= zoom_cutoff), ]

  zoom_corr <- stats::cor(zoom_df$run_a, zoom_df$run_b, method='pearson')
  zoom_corr <- base::round(zoom_corr, digits=3)


  ## The smallest non-zero value from ESPRESSO is 0.01
  ## Convert zero values to 0.001 to show on the log transformed plot
  fake_zero <- 0.001
  plot_df$run_a[plot_df$run_a < 0.01] <- fake_zero
  plot_df$run_b[plot_df$run_b < 0.01] <- fake_zero

  max_value <- base::max(base::max(plot_df$run_a),
                         base::max(plot_df$run_b))

  highest_tick <- base::ceiling(base::log10(max_value))
  break_powers <- -3:highest_tick
  breaks <- 10^break_powers
  last_break <- breaks[base::length(breaks)]
  limits <- base::c(breaks[1], last_break)
  labels <- base::sapply(break_powers, function(x) {
    if (x == -1) {
      return('0.1')
    }
    if (x == 0) {
      return('1')
    }
    if (x == 1) {
      return('10')
    }
    return(base::paste0('10^', x))
  })
  ## 0.001 values are actually zero
  labels[1] <- '0'

  ## Input can be either isoform or gene
  num_isoforms <- base::nrow(plot_df)
  title <- base::paste0(sample, ': abundance (n=', num_isoforms, ', r=', corr, ')')
  plot <- ggplot2::ggplot(data=plot_df) +
    ggplot2::geom_point(ggplot2::aes(x=run_a, y=run_b), color=colors$blue,
                        alpha=0.1) +
    ggplot2::geom_line(ggplot2::aes(x=x, y=x),
                       base::data.frame(x=base::c(fake_zero, max_value)),
                       color=colors$black) +
    ggplot2::scale_x_continuous(transform='log10', breaks=breaks, labels=labels,
                                limits=limits, guide='axis_logticks') +
    ggplot2::scale_y_continuous(transform='log10', breaks=breaks, labels=labels,
                                limits=limits, guide='axis_logticks') +
    create_ggplot_theme() +
    ggplot2::labs(x=run_a_name, y=run_b_name, title=title)
  save_plot(plot, sample, out_dir)

  zoom_breaks <- base::seq(from=0, to=zoom_cutoff, by=100)
  zoom_labels <- base::paste0(zoom_breaks)
  zoom_limits <- base::c(0, zoom_cutoff)
  zoom_num_isoforms <- base::nrow(zoom_df)
  zoom_title <- base::paste0(sample, ': isoform abundance <= ', zoom_cutoff,
                             ' (n=', zoom_num_isoforms, ', r=', zoom_corr, ')')
  plot <- ggplot2::ggplot(data=zoom_df) +
    ggplot2::geom_point(ggplot2::aes(x=run_a, y=run_b), color=colors$blue,
                        alpha=0.1) +
    ggplot2::geom_line(ggplot2::aes(x=x, y=x),
                       base::data.frame(x=base::c(0, zoom_cutoff)),
                       color=colors$black) +
    ggplot2::scale_x_continuous(breaks=zoom_breaks, labels=zoom_labels,
                                limits=zoom_limits) +
    ggplot2::scale_y_continuous(breaks=zoom_breaks, labels=zoom_labels,
                                limits=zoom_limits) +
    create_ggplot_theme() +
    ggplot2::labs(x=run_a_name, y=run_b_name, title=zoom_title)
  save_plot(plot, base::paste0(sample, '_', zoom_cutoff), out_dir)

  return(corr)
}

make_plot_for_all_samples <- function(samples, corr_values, out_dir) {
  colors <- get_colors()
  num_samples <- base::length(samples)
  title <- base::paste0('correlation by sample (n=', num_samples, ')')
  df <- base::data.frame(x=1, y=corr_values)

  min_value <- base::min(corr_values)
  max_value <- base::max(corr_values)
  min_break <- base::floor(min_value * 100) / 100
  max_break <- base::ceiling(max_value * 100) / 100
  breaks <- base::seq(from=min_break, to=max_break, by=0.01)
  limits <- base::c(min_break, max_break)

  plot <- ggplot2::ggplot(data=df) +
    ggplot2::geom_violin(ggplot2::aes(x=x, y=y),
                         color=colors$orange) +
    ggplot2::geom_point(ggplot2::aes(x=x, y=y), color=colors$blue,
                        alpha=0.5,
                        position=ggplot2::position_jitter(width=0.1, height=0)) +
    ggplot2::scale_x_continuous(breaks=NULL) +
    ggplot2::scale_y_continuous(breaks=breaks, limits=limits) +
    create_ggplot_theme() +
    ggplot2::labs(x='samples', y='correlation', title=title)
  save_plot(plot, 'all', out_dir)
}

create_out_dir <- function(path) {
  if (base::dir.exists(path)) {
    return()
  }
  base::cat('mkdir -p ', path, '\n')
  base::dir.create(path, recursive=TRUE)
}

main <- function() {
  create_out_dir(out_dir)
  df_and_sample_names <- read_comparison_tsv(comparison_tsv_path)
  df <- df_and_sample_names$df
  samples <- df_and_sample_names$sample_names
  corr_values <- base::vector(mode='numeric')
  for (sample in samples) {
    corr <- make_plot_for_sample(df, sample, run_a_name, run_b_name, out_dir)
    corr_values <- base::append(corr_values, corr)
  }

  make_plot_for_all_samples(samples, corr_values, out_dir)
}

main()
