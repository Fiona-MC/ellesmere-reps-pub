# plot_composition_combined.r
# Eukaryote (top) and Microbe (bottom) composition errorbars stacked vertically,
# no titles. Shares legend across both panels.
# Output: final_plots/composition_errorbars_combined.pdf
# Run from the ellesmere/ project root.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

genus_colors <- unname(palette.colors(palette = "Polychrome 36"))[-1]

rsa <- read.csv("data/FMC_processed/replicate_sample_assignments_newNames.csv",
                check.names = FALSE)
rsa <- rsa[, c("library_id", "library_type")]
rsa <- rsa[!duplicated(rsa$library_id), ]

loc_depth_num <- function(x) as.numeric(sub("m.*", "", x))

make_panel <- function(prefix, label) {
    ss_path <- file.path("data/FMC_processed", paste0(prefix, "_ss.csv"))
    ds_path <- file.path("data/FMC_processed", paste0(prefix, "_ds.csv"))
    data <- rbind(read.csv(ss_path, na.strings = ""),
                  read.csv(ds_path, na.strings = ""))
    data <- data[!is.na(data$genus), ]
    data <- data[, setdiff(names(data), c("tech_replicate", "library_type"))]
    data <- merge(data, rsa, by = "library_id", all.x = TRUE)
    data <- data[!is.na(data$profile_location) & !is.na(data$core_replicate), ]
    data <- data[data$profile_location != "54m", ]

    locs <- unique(data$profile_location)
    strat_levels <- locs[order(loc_depth_num(locs), locs)]

    data_pct <- data %>%
        group_by(library_id, profile_location) %>%
        mutate(total_reads = sum(n_reads, na.rm = TRUE),
               pct         = ifelse(total_reads > 0, n_reads / total_reads * 100, NA_real_)) %>%
        ungroup()

    stats <- data_pct %>%
        group_by(profile_location, genus) %>%
        summarise(mean_pct = mean(pct, na.rm = TRUE),
                  sd_pct   = sd(pct,  na.rm = TRUE),
                  n_reps   = sum(!is.na(pct)), .groups = "drop") %>%
        mutate(sd_pct = ifelse(is.na(sd_pct), 0, sd_pct),
               ci95   = ifelse(n_reps > 1, 2 * sd_pct / sqrt(n_reps), 0))

    genus_order <- stats %>%
        group_by(genus) %>%
        summarise(tot = sum(mean_pct), .groups = "drop") %>%
        arrange(tot) %>% pull(genus)

    top10 <- tail(genus_order, 10)

    stats <- stats %>%
        mutate(genus            = factor(genus, levels = genus_order),
               profile_location = factor(profile_location, levels = strat_levels, ordered = TRUE)) %>%
        arrange(profile_location, desc(as.integer(genus))) %>%
        group_by(profile_location) %>%
        mutate(cum_top = cumsum(mean_pct), cum_bottom = cum_top - mean_pct) %>%
        ungroup()

    taxa_colors <- rep(genus_colors, length.out = nlevels(stats$genus))
    names(taxa_colors) <- levels(stats$genus)
    eb_data <- stats %>% filter(genus %in% top10)

    ggplot(stats, aes(x = profile_location, fill = genus)) +
        geom_col(aes(y = mean_pct), position = "stack", width = 0.75) +
        geom_errorbar(data = eb_data,
                      aes(ymin = cum_top - ci95, ymax = cum_top),
                      width = 0.3, linewidth = 0.5, color = "grey15",
                      position = "identity") +
        geom_point(data = eb_data, aes(y = cum_top),
                   size = 1.2, color = "grey15", position = "identity") +
        scale_fill_manual(values = taxa_colors, name = label) +
        labs(x = "Profile location", y = "% reads") +
        theme_bw(base_size = 13) +
        theme(axis.text.x    = element_text(angle = 45, hjust = 1),
              legend.position = "right")
}

p_euk <- make_panel("ellesmere_long",     "Eukaryote")
p_mic <- make_panel("ellesmere_mic_long", "Microbe")

combined <- p_euk / p_mic

out_path <- "final_plots/composition_errorbars_combined.pdf"
ggsave(out_path, combined, width = 14, height = 12)
ggsave(sub("\\.pdf$", ".png", out_path), combined, width = 14, height = 12, dpi = 300)
cat("Saved:", out_path, "\n")
