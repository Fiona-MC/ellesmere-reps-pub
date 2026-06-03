# run_aitchison_distances.r
# Computes Aitchison distances between all library pairs, classifies pairs as
# tech-rep / bio-rep / non-rep, and runs permutation tests (9999 perms) for
# all pairwise group comparisons. Saves plot-ready data to RDS so that
# plot_aitchison_distances.r can be re-run without repeating this computation.
#
# Pair classification (within same library_type for ss/ds subsets; any for combined):
#   Tech rep  : same profile_location + same library_type + same core_replicate + diff library_id
#   Core rep  : same profile_location + same library_type + diff core_replicate
#   Non-rep   : diff profile_location
#   Cross-type: same profile_location + diff library_type → included by default (exclude_cross_type = FALSE)
#
# Aitchison distance = Euclidean on CLR-transformed proportions (+0.5 pseudocount).
#
# Output: final_results/plot_data/aitchison_distances.rds
# Run from the ellesmere/ project root.

library(vegan)
library(dplyr)

rsa <- read.csv("data/FMC_processed/replicate_sample_assignments_newNames.csv",
                check.names = FALSE)
rsa <- rsa[, c("library_id", "library_type")]
rsa <- rsa[!duplicated(rsa$library_id), ]

data_configs <- list(
    list(stem = "ellesmere_long",     label = "Eukaryote"),
    list(stem = "ellesmere_mic_long", label = "Microbe")
)

sig_stars <- function(p) {
    if (is.na(p))  return("NA")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    return("ns")
}

classify_pairs <- function(d_mat, lib_meta, lt_filter = NULL, exclude_cross_type = FALSE) {
    if (!is.null(lt_filter)) {
        lib_meta <- lib_meta[lib_meta$library_type == lt_filter, ]
        lib_meta <- lib_meta[lib_meta$library_id %in% rownames(d_mat), ]
    }
    libs <- lib_meta$library_id
    n    <- length(libs)
    if (n < 2) return(data.frame())
    rows <- vector("list", n * (n - 1) / 2)
    k    <- 0L
    for (i in seq_len(n - 1)) {
        for (j in seq(i + 1, n)) {
            li <- libs[i]
            lj <- libs[j]
            if (!li %in% rownames(d_mat) || !lj %in% rownames(d_mat)) next
            mi <- lib_meta[lib_meta$library_id == li, ]
            mj <- lib_meta[lib_meta$library_id == lj, ]
            same_hole <- mi$profile_location == mj$profile_location
            same_lt   <- mi$library_type == mj$library_type
            same_br   <- mi$core_replicate == mj$core_replicate
            pair_type <- if (!same_hole) {
                "Non-rep"
            } else if (same_lt || !exclude_cross_type) {
                if (same_br) "Tech rep" else "Core rep"
            } else {
                NA_character_
            }
            if (is.na(pair_type)) next
            k <- k + 1L
            rows[[k]] <- data.frame(pair_type = pair_type, dist = d_mat[li, lj],
                                    stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, rows[seq_len(k)])
}

all_summary <- list()
all_raw     <- list()

for (cfg in data_configs) {
    cat("\n===", cfg$label, "===\n")
    counts <- rbind(
        read.csv(file.path("data/FMC_processed", paste0(cfg$stem, "_ss.csv")), na.strings = ""),
        read.csv(file.path("data/FMC_processed", paste0(cfg$stem, "_ds.csv")), na.strings = "")
    )
    counts <- counts[!is.na(counts$genus), ]
    counts <- counts[, setdiff(names(counts), c("tech_replicate", "library_type"))]
    counts <- merge(counts, rsa, by = "library_id", all.x = TRUE)
    counts <- counts[!is.na(counts$profile_location), ]
    counts <- counts[counts$profile_location != "54m", ]

    lib_wide <- counts %>%
        group_by(library_id, genus) %>%
        summarise(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = genus, values_from = n_reads, values_fill = 0)

    mat <- as.matrix(lib_wide[, -1])
    rownames(mat) <- lib_wide$library_id
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]

    pseudo  <- mat + 0.5
    prop_ps <- pseudo / rowSums(pseudo)
    clr_mat <- log(prop_ps) - rowMeans(log(prop_ps))
    d_mat   <- as.matrix(dist(clr_mat, method = "euclidean"))

    lib_meta <- counts[!duplicated(counts$library_id),
                       c("library_id", "library_type", "profile_location", "core_replicate")]
    lib_meta <- lib_meta[lib_meta$library_id %in% rownames(d_mat), ]

    for (subset in c("combined", "ss", "ds")) {
        lt_filter <- if (subset == "combined") NULL else subset
        pairs     <- classify_pairs(d_mat, lib_meta, lt_filter)
        if (nrow(pairs) == 0) next
        summ <- pairs %>%
            group_by(pair_type) %>%
            summarise(mean_dist = mean(dist, na.rm = TRUE),
                      se_dist   = sd(dist, na.rm = TRUE) / sqrt(sum(!is.na(dist))),
                      n         = sum(!is.na(dist)), .groups = "drop") %>%
            mutate(dataset = cfg$label, subset = subset)
        all_summary[[length(all_summary) + 1]] <- summ
        all_raw[[paste(cfg$label, subset, sep = ":::")]] <-
            mutate(pairs, dataset = cfg$label, subset = subset)
    }
}

plot_df <- do.call(rbind, all_summary)
plot_df$pair_type <- factor(plot_df$pair_type, levels = c("Tech rep", "Core rep", "Non-rep"))
plot_df$subset    <- factor(plot_df$subset,
                            levels = c("combined", "ss", "ds"),
                            labels = c("ss + ds combined", "ss only", "ds only"))
plot_df$dataset   <- factor(plot_df$dataset, levels = c("Eukaryote", "Microbe"))

# ── Permutation tests for all pairwise comparisons ────────────────────────────

subset_labels <- c("combined" = "ss + ds combined", "ss" = "ss only", "ds" = "ds only")
dodge_w  <- 0.7
n_levels <- 3
offsets  <- (seq_len(n_levels) - (n_levels + 1) / 2) * dodge_w / n_levels
names(offsets) <- c("Tech rep", "Core rep", "Non-rep")
x_pos <- c(Eukaryote = 1, Microbe = 2)

comparisons  <- list(c("Tech rep", "Core rep"), c("Core rep", "Non-rep"), c("Tech rep", "Non-rep"))
panel_top    <- max(plot_df$mean_dist + plot_df$se_dist, na.rm = TRUE)
bracket_step <- panel_top * 0.07

sig_rows <- list()
for (ds_label in c("Eukaryote", "Microbe")) {
    for (sub in c("combined", "ss", "ds")) {
        key <- paste(ds_label, sub, sep = ":::")
        raw <- all_raw[[key]]
        if (is.null(raw)) next
        top_y <- plot_df %>%
            filter(dataset == ds_label, subset == subset_labels[sub]) %>%
            summarise(top = max(mean_dist + se_dist, na.rm = TRUE)) %>%
            pull(top)
        for (k in seq_along(comparisons)) {
            pt_a <- comparisons[[k]][1]; pt_b <- comparisons[[k]][2]
            vals_a <- raw$dist[raw$pair_type == pt_a]
            vals_b <- raw$dist[raw$pair_type == pt_b]
            if (length(vals_a) < 2 || length(vals_b) < 2) next
            obs_diff <- abs(mean(vals_a) - mean(vals_b))
            pooled   <- c(vals_a, vals_b)
            n_a      <- length(vals_a)
            perm_diffs <- replicate(9999, {
                s <- sample(pooled)
                abs(mean(s[seq_len(n_a)]) - mean(s[seq(n_a + 1, length(s))]))
            })
            p_val <- (sum(perm_diffs >= obs_diff) + 1) / (9999 + 1)
            x_ctr <- x_pos[ds_label]
            x_l   <- x_ctr + offsets[pt_a]
            x_r   <- x_ctr + offsets[pt_b]
            y_br  <- top_y + bracket_step * k
            sig_rows[[length(sig_rows) + 1]] <- data.frame(
                subset  = subset_labels[sub], dataset = ds_label,
                x_left  = min(x_l, x_r), x_right = max(x_l, x_r),
                x_mid   = (x_l + x_r) / 2,
                y_bar   = y_br, y_tick = y_br - bracket_step * 0.3,
                label   = sig_stars(p_val), stringsAsFactors = FALSE)
        }
    }
}

sig_df         <- do.call(rbind, sig_rows)
sig_df$subset  <- factor(sig_df$subset, levels = c("ss + ds combined", "ss only", "ds only"))
sig_df$dataset <- factor(sig_df$dataset, levels = c("Eukaryote", "Microbe"))

# ── Save ──────────────────────────────────────────────────────────────────────

out_dir <- "final_results/plot_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(list(plot_df = plot_df, sig_df = sig_df, raw = all_raw),
        file.path(out_dir, "aitchison_distances.rds"))
cat("Saved: final_results/plot_data/aitchison_distances.rds\n")

write.csv(
    plot_df %>% select(dataset, subset, pair_type, mean_dist, se_dist, n),
    "final_results/aitchison_summary.csv", row.names = FALSE
)
cat("Saved: final_results/aitchison_summary.csv\n")
cat("\nDone.\n")
