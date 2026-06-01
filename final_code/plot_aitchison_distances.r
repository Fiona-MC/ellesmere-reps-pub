# plot_aitchison_distances.r
# Bar plot of mean Aitchison distance for tech-rep, bio-rep, and non-rep pairs.
# Reads pre-computed plot data saved by run_aitchison_distances.r.
#
# Requires: final_results/plot_data/aitchison_distances.rds
#   (run run_aitchison_distances.r first)
#
# Output: final_plots/aitchison_distances.pdf
# Run from the ellesmere/ project root.

library(ggplot2)

rds_path <- "final_results/plot_data/aitchison_distances.rds"
if (!file.exists(rds_path)) stop("Missing: ", rds_path, " — run run_aitchison_distances.r first")

dat     <- readRDS(rds_path)
plot_df <- dat$plot_df
sig_df  <- dat$sig_df

pair_colors <- c("Tech rep" = "#59A14F", "Core rep" = "#F28E2B", "Non-rep" = "#4E79A7")
dodge_w     <- 0.7

out_dir <- "final_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

p <- ggplot(plot_df, aes(x = dataset, y = mean_dist, fill = pair_type)) +
    geom_col(position = position_dodge(width = dodge_w), width = 0.6) +
    geom_errorbar(aes(ymin = mean_dist - se_dist, ymax = mean_dist + se_dist),
                  position = position_dodge(width = dodge_w),
                  width = 0.2, linewidth = 0.5) +
    geom_segment(data = sig_df,
                 aes(x = x_left, xend = x_right, y = y_bar, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.45, colour = "grey20") +
    geom_segment(data = sig_df,
                 aes(x = x_left, xend = x_left, y = y_tick, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.45, colour = "grey20") +
    geom_segment(data = sig_df,
                 aes(x = x_right, xend = x_right, y = y_tick, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.45, colour = "grey20") +
    geom_text(data = sig_df,
              aes(x = x_mid, y = y_bar + (y_bar - y_tick) * 0.5, label = label),
              inherit.aes = FALSE, size = 3, colour = "grey10") +
    scale_fill_manual(values = pair_colors, name = "Pair type",
                      labels = c("Tech rep" = "Tech rep (same core)",
                                 "Core rep" = "Core rep (same location, diff core)",
                                 "Non-rep"  = "Non-rep (diff location)")) +
    facet_wrap(~ subset, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x       = NULL,
         y       = expression("Mean Aitchison distance (" %+-% "SE)"),
         title   = "Replicate vs non-replicate distances  —  Aitchison",
         caption = expression("Brackets: permutation test (9999 perms).   *** p<0.001,  ** p<0.01,  * p<0.05,  ns p" >= 0.05)) +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          strip.text       = element_text(size = 11),
          legend.key.size  = unit(0.45, "cm"),
          plot.caption     = element_text(size = 8, colour = "grey40"))

ggsave(file.path(out_dir, "aitchison_distances.pdf"), p, width = 10, height = 5)
ggsave(file.path(out_dir, "aitchison_distances.png"), p, width = 10, height = 5, dpi = 300)
cat("Saved: final_plots/aitchison_distances.pdf/.png\n")
