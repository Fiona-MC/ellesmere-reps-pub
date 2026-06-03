# run_ss_vs_ds.r
# Tests whether ss and ds libraries differ in composition after controlling for profile_location.
#
# Approach: within-location matched pairs using Aitchison distance.
#   Reads are aggregated to (profile_location, library_type, core_replicate) before computing distances,
#   so each core replicate contributes one profile per library type.
#   Only locations with ≥ 2 core_replicate labels present in BOTH ss and ds are used.
#
#   Three pair types compared:
#     ss-ss         : ssA vs ssB  (core-rep variation within ss)
#     ds-ds         : dsA vs dsB  (core-rep variation within ds)
#     matched-cross : ssX vs dsX  (library-type distance, same core-rep)
#
#   If matched-cross ≈ ss-ss and ds-ds, ss and ds are not more different from each
#   other than core replicates within the same library type.
#
#   Statistical test: lmer(distance ~ pair_type + (1|profile_location)), likelihood ratio ANOVA.
#   Pairwise contrasts computed from REML model.
#
# Outputs:
#   final_results/ss_vs_ds_contrasts.csv   — pairwise contrast table
#   final_results/ss_vs_ds_means.csv       — mean distances per group
# Run from the ellesmere/ project root.

library(lme4)
library(lmerTest)
library(dplyr)

rsa <- read.csv("data/FMC_processed/replicate_sample_assignments_newNames.csv",
                check.names = FALSE)
rsa <- rsa[, c("library_id", "library_type")]
rsa <- rsa[!duplicated(rsa$library_id), ]

data_configs <- list(
    list(stem = "ellesmere_long",     label = "Euk"),
    list(stem = "ellesmere_mic_long", label = "Mic")
)

out_dir <- "final_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

all_contrasts <- list()
all_means     <- list()
all_pairs     <- list()

for (cfg in data_configs) {
    cat("\n========================================\n", cfg$label, "\n")

    counts <- rbind(
        read.csv(file.path("data/FMC_processed", paste0(cfg$stem, "_ss.csv")), na.strings = ""),
        read.csv(file.path("data/FMC_processed", paste0(cfg$stem, "_ds.csv")), na.strings = "")
    )
    counts <- counts[!is.na(counts$genus), ]
    counts <- counts[, setdiff(names(counts), c("tech_replicate", "library_type"))]
    counts <- merge(counts, rsa, by = "library_id", all.x = TRUE)
    counts <- counts[!is.na(counts$profile_location) & !is.na(counts$core_replicate), ]
    counts <- counts[counts$profile_location != "54m", ]

    # Aggregate reads to (profile_location, library_type, core_replicate)
    agg <- counts %>%
        group_by(profile_location, library_type, core_replicate, genus) %>%
        summarise(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop") %>%
        mutate(sample_id = paste(profile_location, library_type, core_replicate, sep = ":::"))

    wide <- agg %>%
        select(sample_id, genus, n_reads) %>%
        tidyr::pivot_wider(names_from = genus, values_from = n_reads, values_fill = 0)

    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$sample_id
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]

    # Aitchison distance: CLR on counts + 0.5 pseudocount, then Euclidean
    pseudo  <- mat + 0.5
    prop_ps <- pseudo / rowSums(pseudo)
    clr_mat <- log(prop_ps) - rowMeans(log(prop_ps))
    d_mat   <- as.matrix(dist(clr_mat, method = "euclidean"))

    samp_meta <- data.frame(
        sample_id        = rownames(mat),
        profile_location = sub(":::.*", "", rownames(mat)),
        library_type     = sub(".*:::(ss|ds):::.*", "\\1", rownames(mat)),
        core_rep         = sub(".*:::(ss|ds):::(.*)", "\\2", rownames(mat)),
        stringsAsFactors = FALSE
    )

    # Build pairs: only locations with ≥ 2 matched core-reps in both ss and ds
    pairs_list <- list()
    for (h in unique(samp_meta$profile_location)) {
        h_meta <- samp_meta[samp_meta$profile_location == h, ]
        matched_brs <- intersect(
            h_meta$core_rep[h_meta$library_type == "ss"],
            h_meta$core_rep[h_meta$library_type == "ds"]
        )
        if (length(matched_brs) < 2) next

        h_sub  <- h_meta[h_meta$core_rep %in% matched_brs, ]
        get_id <- function(lt, br)
            h_sub$sample_id[h_sub$library_type == lt & h_sub$core_rep == br]

        # matched-cross: ssX vs dsX for each matched core-rep
        for (br in matched_brs) {
            id_ss <- get_id("ss", br); id_ds <- get_id("ds", br)
            if (length(id_ss) == 1 && length(id_ds) == 1 &&
                id_ss %in% rownames(d_mat) && id_ds %in% rownames(d_mat)) {
                pairs_list[[length(pairs_list) + 1]] <- data.frame(
                    profile_location = h, pair_type = "matched-cross", dist = d_mat[id_ss, id_ds])
            }
        }

        # within-type bio-rep pairs
        br_combos <- combn(matched_brs, 2, simplify = FALSE)
        for (brs in br_combos) {
            for (lt in c("ss", "ds")) {
                id1 <- get_id(lt, brs[1]); id2 <- get_id(lt, brs[2])
                if (length(id1) == 1 && length(id2) == 1 &&
                    id1 %in% rownames(d_mat) && id2 %in% rownames(d_mat)) {
                    pairs_list[[length(pairs_list) + 1]] <- data.frame(
                        profile_location = h, pair_type = paste0(lt, "-", lt), dist = d_mat[id1, id2])
                }
            }
        }
    }

    pairs <- do.call(rbind, pairs_list)

    # Keep only holes that have all three pair types
    valid_holes <- Reduce(intersect, list(
        unique(pairs$profile_location[pairs$pair_type == "ss-ss"]),
        unique(pairs$profile_location[pairs$pair_type == "ds-ds"]),
        unique(pairs$profile_location[pairs$pair_type == "matched-cross"])
    ))
    pairs <- pairs[pairs$profile_location %in% valid_holes, ]

    cat("  Locations used:", length(valid_holes), ":", paste(sort(valid_holes), collapse = ", "), "\n")
    cat("Pair counts:\n"); print(table(pairs$pair_type))

    pairs$pair_type <- factor(pairs$pair_type,
                              levels = c("ss-ss", "matched-cross", "ds-ds"))

    means_df <- pairs %>%
        group_by(pair_type) %>%
        summarise(mean_dist = mean(dist, na.rm = TRUE),
                  se_dist   = sd(dist, na.rm = TRUE) / sqrt(n()),
                  n         = n(), .groups = "drop") %>%
        mutate(dataset = cfg$label)
    cat("\nMean Aitchison distance by pair type:\n"); print(means_df)
    all_means[[cfg$label]] <- means_df

    # Likelihood ratio test
    model_full    <- lmer(dist ~ pair_type + (1 | profile_location), data = pairs, REML = FALSE)
    model_null    <- lmer(dist ~ 1            + (1 | profile_location), data = pairs, REML = FALSE)
    lrt           <- anova(model_null, model_full)
    cat("\nLikelihood ratio test (pair_type effect):\n"); print(lrt)

    # Pairwise contrasts from REML model
    model_reml <- lmer(dist ~ pair_type + (1 | profile_location), data = pairs, REML = TRUE)
    contr_mat  <- rbind(
        "matched-cross vs ss-ss" = c(0,  1,  0),
        "ds-ds vs ss-ss"         = c(0,  0,  1),
        "matched-cross vs ds-ds" = c(0,  1, -1)
    )
    est    <- contr_mat %*% fixef(model_reml)
    se_c   <- sqrt(diag(contr_mat %*% as.matrix(vcov(model_reml)) %*% t(contr_mat)))
    tval   <- est / se_c
    coef_d <- summary(model_reml)$coefficients
    avg_df <- if (nrow(coef_d) >= 3) mean(coef_d[2:3, "df"]) else coef_d[2, "df"]
    pval   <- 2 * pt(abs(tval), df = avg_df, lower.tail = FALSE)

    contrasts_df <- data.frame(
        dataset  = cfg$label,
        contrast = rownames(est),
        estimate = round(est,  4),
        se       = round(se_c, 4),
        t        = round(tval, 3),
        p        = round(pval, 4)
    )
    cat("\nPairwise contrasts:\n"); print(contrasts_df)
    all_contrasts[[cfg$label]] <- contrasts_df
    all_pairs[[cfg$label]]     <- pairs
}

contrasts_out <- do.call(rbind, all_contrasts)
means_out     <- do.call(rbind, all_means)

write.csv(contrasts_out, file.path(out_dir, "ss_vs_ds_contrasts.csv"), row.names = FALSE)
write.csv(means_out,     file.path(out_dir, "ss_vs_ds_means.csv"),     row.names = FALSE)
cat("\nSaved: final_results/ss_vs_ds_contrasts.csv\n")
cat("Saved: final_results/ss_vs_ds_means.csv\n")

# Save raw pairs and model summaries so plots can be made later without re-running
saveRDS(list(pairs     = all_pairs,
             contrasts = contrasts_out,
             means     = means_out),
        file.path(out_dir, "plot_data/ss_vs_ds.rds"))
cat("Saved: final_results/plot_data/ss_vs_ds.rds\n")
cat("\nDone.\n")
