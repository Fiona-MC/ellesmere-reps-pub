# plot_varpart_barplot.r
# Horizontal stacked bar chart, one bar per taxon, showing the proportion of
# variance explained by each component: ds vs ss / Profile location / Core rep / Residual.
# Taxa sorted by Profile location variance (ascending) so spatially structured taxa are at top.
# One plot per organism type (Euk, Mic).
#
# Requires: final_results/plot_data/var_partitioning_{Euk,Mic}.rds
#   (run run_variance_partitioning.r first)
#
# Outputs:
#   final_plots/varpart_barplot_euk.pdf
#   final_plots/varpart_barplot_mic.pdf
# Run from the ellesmere/ project root.

library(ggplot2)
library(dplyr)
library(tidyr)

comp_colors <- c(`ds vs ss`          = "#984EA3",
                 `Profile location`  = "#4E79A7",
                 `Core rep`          = "#F28E2B",
                 Residual            = "#59A14F")
comp_levels <- c("ds vs ss", "Profile location", "Core rep", "Residual")

data_configs <- list(
    list(label = "Eukaryote", rds = "final_results/plot_data/var_partitioning_Euk.rds"),
    list(label = "Microbe",   rds = "final_results/plot_data/var_partitioning_Mic.rds")
)

out_dir <- "final_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (cfg in data_configs) {
    if (!file.exists(cfg$rds)) {
        stop("Missing: ", cfg$rds, " — run run_variance_partitioning.r first")
    }

    vp     <- readRDS(cfg$rds)
    var_df <- vp$var_df

    # Drop rows where all components are NA (model didn't converge)
    complete_rows <- rowSums(!is.na(var_df[, c("profile_location", "library_type")])) > 0
    var_df <- var_df[complete_rows, ]

    var_long <- var_df %>%
        pivot_longer(cols = c(library_type, profile_location, loc_core_rep, Residual),
                     names_to  = "component",
                     values_to = "variance") %>%
        filter(!is.na(variance)) %>%
        mutate(component = factor(component,
            levels = c("library_type", "profile_location", "loc_core_rep", "Residual"),
            labels = comp_levels))

    # Sort by profile_location variance (ascending) so high-location taxa appear at top of plot
    genus_order <- var_df %>%
        mutate(hole_var = ifelse(is.na(profile_location), 0, profile_location)) %>%
        arrange(hole_var) %>%
        pull(genus)

    var_long <- var_long %>%
        mutate(genus = factor(genus, levels = genus_order))

    n_genera   <- length(genus_order)
    bar_height <- max(0.4, min(0.85, 20 / n_genera))
    fig_height <- max(4, n_genera * 0.22 + 1.5)

    p <- ggplot(var_long, aes(x = variance, y = genus, fill = component)) +
        geom_col(position = "stack", width = bar_height) +
        scale_fill_manual(values = comp_colors, name = "Component",
                          breaks = comp_levels) +
        scale_x_continuous(labels = scales::percent, limits = c(0, 1),
                           expand = c(0, 0)) +
        labs(x     = "Proportion of variance",
             y     = NULL,
             title = paste0(cfg$label, " — variance partitioning per taxon")) +
        theme_bw(base_size = 10) +
        theme(axis.text.y          = element_text(size = max(5, min(9, 160 / n_genera))),
              legend.position      = "bottom",
              legend.key.size      = unit(0.4, "cm"),
              panel.grid.major.y   = element_blank())

    out_path <- file.path(out_dir,
                          paste0("varpart_barplot_", tolower(cfg$label), ".pdf"))
    ggsave(out_path, p, width = 8, height = fig_height)
    ggsave(sub("\\.pdf$", ".png", out_path), p, width = 8, height = fig_height, dpi = 300)
    cat("Saved:", out_path, "\n")
}

cat("\nDone.\n")
