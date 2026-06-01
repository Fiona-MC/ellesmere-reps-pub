# plot_varpart_violin.r
# Violin + box plots of per-taxon variance components from the combined ss+ds model.
# Components: ds vs ss (library_type fixed effect) / Profile location / Core rep / Residual.
# Faceted by Euk and Mic.
#
# Requires: final_results/plot_data/var_partitioning_{Euk,Mic}.rds
#   (run run_variance_partitioning.r first)
#
# Output: final_plots/varpart_violin.pdf
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

all_var <- list()
for (cfg in data_configs) {
    if (!file.exists(cfg$rds)) {
        stop("Missing: ", cfg$rds, " — run run_variance_partitioning.r first")
    }
    vp     <- readRDS(cfg$rds)
    var_df <- vp$var_df

    var_long <- var_df %>%
        select(genus, library_type, profile_location, loc_core_rep, Residual) %>%
        pivot_longer(cols = c(library_type, profile_location, loc_core_rep, Residual),
                     names_to  = "component",
                     values_to = "variance") %>%
        filter(!is.na(variance)) %>%
        mutate(
            component = factor(component,
                levels = c("library_type", "profile_location", "loc_core_rep", "Residual"),
                labels = comp_levels),
            dataset = cfg$label
        )
    all_var[[cfg$label]] <- var_long
}

combined        <- do.call(rbind, all_var)
combined$dataset <- factor(combined$dataset, levels = c("Eukaryote", "Microbe"))

out_dir <- "final_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

p <- ggplot(combined, aes(x = component, y = variance, fill = component)) +
    geom_violin(alpha = 0.65, trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.8) +
    scale_fill_manual(values = comp_colors, breaks = comp_levels) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ dataset, nrow = 1) +
    labs(x        = NULL,
         y        = "Proportion of variance",
         title    = "Variance partitioning — ss + ds combined",
         subtitle = "Model: read_proportions ~ library_type + (1|profile_location) + (1|profile_location:core_rep)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          strip.text      = element_text(size = 12))

ggsave(file.path(out_dir, "varpart_violin.pdf"), p, width = 10, height = 5)
ggsave(file.path(out_dir, "varpart_violin.png"), p, width = 10, height = 5, dpi = 300)
cat("Saved: final_plots/varpart_violin.pdf/.png\n")
