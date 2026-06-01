# plot_varpart_combined.r
# Eukaryote (left) and Microbe (right) varpart barplots side by side, no titles.
# Shares legend across both panels.
# Output: final_plots/varpart_barplot_combined.pdf
# Run from the ellesmere/ project root.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

comp_colors <- c(`ds vs ss`          = "#984EA3",
                 `Profile location`  = "#4E79A7",
                 `Core rep`          = "#F28E2B",
                 Residual            = "#59A14F")
comp_levels <- c("ds vs ss", "Profile location", "Core rep", "Residual")

data_configs <- list(
    list(label = "Eukaryote", rds = "final_results/plot_data/var_partitioning_Euk.rds"),
    list(label = "Microbe",   rds = "final_results/plot_data/var_partitioning_Mic.rds")
)

make_panel <- function(cfg, show_legend) {
    if (!file.exists(cfg$rds)) stop("Missing: ", cfg$rds, " — run run_variance_partitioning.r first")

    vp     <- readRDS(cfg$rds)
    var_df <- vp$var_df

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

    genus_order <- var_df %>%
        mutate(hole_var = ifelse(is.na(profile_location), 0, profile_location)) %>%
        arrange(hole_var) %>%
        pull(genus)

    var_long <- var_long %>%
        mutate(genus = factor(genus, levels = genus_order))

    n_genera   <- length(genus_order)
    bar_height <- max(0.4, min(0.85, 20 / n_genera))

    p <- ggplot(var_long, aes(x = variance, y = genus, fill = component)) +
        geom_col(position = "stack", width = bar_height) +
        scale_fill_manual(values = comp_colors, name = "Component",
                          breaks = comp_levels) +
        scale_x_continuous(labels = scales::percent, limits = c(0, 1),
                           expand = c(0, 0)) +
        labs(x = "Proportion of variance", y = NULL) +
        theme_bw(base_size = 10) +
        theme(axis.text.y        = element_text(size = max(5, min(9, 160 / n_genera))),
              legend.position    = if (show_legend) "bottom" else "none",
              legend.key.size    = unit(0.4, "cm"),
              panel.grid.major.y = element_blank())
    p
}

p_euk <- make_panel(data_configs[[1]], show_legend = FALSE)
p_mic <- make_panel(data_configs[[2]], show_legend = TRUE)

combined <- p_euk | p_mic

out_dir <- "final_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_path <- file.path(out_dir, "varpart_barplot_combined.pdf")

vp_euk <- readRDS(data_configs[[1]]$rds)
vp_mic <- readRDS(data_configs[[2]]$rds)
n_euk  <- sum(rowSums(!is.na(vp_euk$var_df[, c("profile_location", "library_type")])) > 0)
n_mic  <- sum(rowSums(!is.na(vp_mic$var_df[, c("profile_location", "library_type")])) > 0)
fig_height <- max(4, max(n_euk, n_mic) * 0.22 + 2)

ggsave(out_path, combined, width = 16, height = fig_height)
ggsave(sub("\\.pdf$", ".png", out_path), combined, width = 16, height = fig_height, dpi = 300)
cat("Saved:", out_path, "\n")
