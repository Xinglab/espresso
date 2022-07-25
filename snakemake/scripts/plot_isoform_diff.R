library('ggplot2')
library('scales')

args <- base::commandArgs(trailingOnly=TRUE)
in_tsv_path <- args[1]
out_plot_base_path <- args[2]
condition <- args[3]

create_ggplot_theme <- function() {
    return(
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(color='black'),
                       legend.title=ggplot2::element_blank()
                       )
    )
}

read_input_table <- function(tsv_path) {
    return(utils::read.table(tsv_path, sep='\t', header=TRUE))
}

make_plot <- function(plot_df, out_plot_base_path, condition, is_transcript) {
    plot_point_size <- 0.5
    plot_point_alpha <- 0.5
    plot_width <- 4.5
    plot_height<- 3.5

    df_below_2 <- plot_df[plot_df$abs_diff < 2, ]
    df_at_least_2 <- plot_df[plot_df$abs_diff >= 2, ]
    df_over_100 <- plot_df[plot_df$abs_diff > 100, ]
    df_2_to_100 <- df_at_least_2[df_at_least_2$abs_diff <= 100, ]
    n_below_2 <- base::nrow(df_below_2)
    n_at_least_2 <- base::nrow(df_at_least_2)
    n_over_100 <- base::nrow(df_over_100)
    n_2_to_100 <- base::nrow(df_2_to_100)

    histo_x_lab <- 'max for transcript: abs(read_count_diff)'
    histo_y_lab <- 'count transcripts'
    point_x_lab <- 'nothing'
    point_y_lab <- 'nothing'
    if (!is_transcript) {
        histo_x_lab <- 'abs(read_count_diff)'
        histo_y_lab <- 'count (transcript, sample)'
        point_x_lab <- 'smaller_read_count'
        point_y_lab <- 'larger_read_count'
    }

    ## histogram 0 -> 2 by 0.1
    title <- base::paste0('Abundance: ', condition)
    subtitle <- base::paste0('n(<2)=', n_below_2, ', n(>=2)=', n_at_least_2)
    plot <- ggplot2::ggplot(data=df_below_2,
                            ggplot2::aes(x=abs_diff)) +
        ggplot2::geom_histogram(breaks=base::seq(from=0, to=2, by=0.1), fill='blue') +
        ggplot2::labs(x=histo_x_lab,
                      y=histo_y_lab,
                      title=title, subtitle=subtitle) +
        create_ggplot_theme() +
        ggplot2::guides(x=ggplot2::guide_axis(angle=45)) +
        ggplot2::scale_x_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA)) +
        ggplot2::scale_y_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA))

    out_plot_path <- base::paste0(out_plot_base_path, '_hist_0_to_2.png')
    ggplot2::ggsave(plot=plot, out_plot_path, width=plot_width,
                    height=plot_height)

    ## histogram 2 -> 100 by 1
    title <- base::paste0('Abundance: ', condition)
    subtitle <- base::paste0('n(<2)=', n_below_2, ', n(2,100)=', n_2_to_100,  ', n(>100)=', n_over_100)
    plot <- ggplot2::ggplot(data=df_2_to_100,
                            ggplot2::aes(x=abs_diff)) +
        ggplot2::geom_histogram(breaks=base::seq(from=2, to=100, by=1), fill='blue') +
        ggplot2::labs(x=histo_x_lab,
                      y=histo_y_lab,
                      title=title, subtitle=subtitle) +
        create_ggplot_theme() +
        ggplot2::guides(x=ggplot2::guide_axis(angle=45)) +
        ggplot2::scale_x_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA)) +
        ggplot2::scale_y_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA))

    out_plot_path <- base::paste0(out_plot_base_path, '_hist_2_to_100.png')
    ggplot2::ggsave(plot=plot, out_plot_path, width=plot_width,
                    height=plot_height)

    df_at_least_10 <- plot_df[plot_df$abs_diff >= 10, ]
    df_at_least_10_b_below_1000 <- df_at_least_10[df_at_least_10$larger < 1000, ]
    df_at_least_10_b_at_least_1000 <- df_at_least_10[df_at_least_10$larger >= 1000, ]
    ## (a, b) scatterplot: diff > 10, b < 1000
    title <- base::paste0('Abundance: ', condition)
    subtitle <- base::paste0('n=', base::nrow(df_at_least_10_b_below_1000))
    plot <- ggplot2::ggplot(data=df_at_least_10_b_below_1000,
                            ggplot2::aes(x=smaller, y=larger)) +
        ggplot2::geom_point(color='blue', size=plot_point_size, alpha=plot_point_alpha) +
        ggplot2::labs(x=point_x_lab,
                      y=point_y_lab,
                      title=title, subtitle=subtitle) +
        create_ggplot_theme() +
        ggplot2::guides(x=ggplot2::guide_axis(angle=45)) +
        ggplot2::scale_x_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA)) +
        ggplot2::scale_y_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA))

    out_plot_path <- base::paste0(out_plot_base_path, '_a_b_diff_10_b_below_1000.png')
    ggplot2::ggsave(plot=plot, out_plot_path, width=plot_width,
                    height=plot_height)

    ## (a, b) scatterplot: diff > 10, b >= 1000
    title <- base::paste0('Abundance: ', condition)
    subtitle <- base::paste0('n=', base::nrow(df_at_least_10_b_at_least_1000))
    plot <- ggplot2::ggplot(data=df_at_least_10_b_at_least_1000,
                            ggplot2::aes(x=smaller, y=larger)) +
        ggplot2::geom_point(color='blue', size=plot_point_size, alpha=plot_point_alpha) +
        ggplot2::labs(x=point_x_lab,
                      y=point_y_lab,
                      title=title, subtitle=subtitle) +
        create_ggplot_theme() +
        ggplot2::guides(x=ggplot2::guide_axis(angle=45)) +
        ggplot2::scale_x_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA)) +
        ggplot2::scale_y_continuous(labels=scales::label_number(big.mark=','),
                                    limits=c(0, NA))

    out_plot_path <- base::paste0(out_plot_base_path, '_a_b_diff_10_b_at_least_1000.png')
    ggplot2::ggsave(plot=plot, out_plot_path, width=plot_width,
                    height=plot_height)
}

plot_df <- read_input_table(in_tsv_path)
transcript_df <- plot_df[plot_df$type == 'transcript', ]
sample_df <- plot_df[plot_df$type == 'sample', ]
make_plot(transcript_df, base::paste0(out_plot_base_path, '_transcript'), condition, TRUE)
make_plot(sample_df, base::paste0(out_plot_base_path, '_sample'), condition, FALSE)
