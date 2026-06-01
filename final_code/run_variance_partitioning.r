# run_variance_partitioning.r
# Per-taxon variance partitioning for Ellesmere data, ss and ds combined.
#
# Model: read_proportions ~ library_type + (1|profile_location) + (1|profile_location:core_replicate)
#   library_type (ds vs ss) is a fixed effect; variance = var(X*beta) across observations.
#   (1|profile_location)              captures variance between locations (depth/position).
#   (1|profile_location:core_rep)     captures variance between core replicates within a location.
#   Residual             captures within-bio-rep (tech-rep / lane) variance.
#
# Output: final_results/plot_data/var_partitioning_{Euk,Mic}.rds
# Run from the ellesmere/ project root.

library(lme4)
library(dplyr)

# proportionsCPM: divide each library (column) by its total reads / 1e6
proportionsCPM <- function(count_matrix) {
    col_sums <- colSums(count_matrix, na.rm = TRUE) / 1e6
    sweep(count_matrix, 2, col_sums, "/")
}

# ── Parameters ────────────────────────────────────────────────────────────────
minReadsPerTaxon <- 1000
minSamples       <- 5

data_configs <- list(
    list(prefix = "ellesmere_long",     label = "Euk"),
    list(prefix = "ellesmere_mic_long", label = "Mic")
)

# ── Load replicate structure ──────────────────────────────────────────────────

rsa <- read.csv("data/FMC_processed/replicate_sample_assignments_newNames.csv",
                check.names = FALSE)
rsa <- rsa[, c("library_id", "library_type", "profile_location", "core_replicate")]
rsa <- rsa[!duplicated(rsa$library_id), ]

# ── Core variance partitioning function ──────────────────────────────────────

get_var_explained <- function(reads_table, lib_meta) {

    reads_norm <- proportionsCPM(reads_table)
    n_genera   <- nrow(reads_norm)

    results <- data.frame(
        genus            = rownames(reads_norm),
        library_type     = NA_real_,
        profile_location = NA_real_,
        loc_core_rep     = NA_real_,
        Residual         = NA_real_,
        stringsAsFactors = FALSE
    )

    get_vc <- function(vc_df, grp) {
        v <- vc_df$vcov[vc_df$grp == grp]
        if (length(v) == 0 || is.na(v[1])) 0 else v[1]
    }

    for (i in seq_len(n_genera)) {
        row_vals <- as.numeric(reads_norm[i, ])
        if (var(row_vals, na.rm = TRUE) == 0) {
            results[i, c("library_type", "profile_location", "loc_core_rep", "Residual")] <- 0
            next
        }

        gdata <- data.frame(
            read_proportions       = row_vals,
            library_type     = lib_meta$library_type,
            profile_location = lib_meta$profile_location,
            core_rep         = lib_meta$core_replicate
        )
        gdata$loc_core_rep <- paste(gdata$profile_location, gdata$core_rep, sep = ":")

        has_both_lt    <- length(unique(gdata$library_type)) == 2
        has_multi_loc  <- length(unique(gdata$profile_location)) >= 2
        has_nested     <- any(tapply(gdata$core_rep, gdata$profile_location,
                                     function(x) length(unique(x))) >= 2, na.rm = TRUE)

        tryCatch({
            if (has_both_lt && has_multi_loc && has_nested) {
                model     <- lmer(read_proportions ~ library_type + (1|profile_location) + (1|loc_core_rep),
                                  data = gdata, REML = TRUE)
                vc        <- as.data.frame(VarCorr(model))
                X_beta    <- model.matrix(model) %*% fixef(model)
                var_fixed <- var(as.numeric(X_beta))
                var_loc   <- get_vc(vc, "profile_location")
                var_hbr   <- get_vc(vc, "loc_core_rep")
                var_resid <- get_vc(vc, "Residual")
                total_var <- var_fixed + var_loc + var_hbr + var_resid
                if (total_var == 0) next
                results[i, "library_type"]     <- var_fixed / total_var
                results[i, "profile_location"] <- var_loc   / total_var
                results[i, "loc_core_rep"]     <- var_hbr   / total_var
                results[i, "Residual"]         <- var_resid / total_var

            } else if (has_multi_loc) {
                model     <- lmer(read_proportions ~ (1|profile_location), data = gdata, REML = TRUE)
                vc        <- as.data.frame(VarCorr(model))
                var_loc   <- get_vc(vc, "profile_location")
                var_resid <- get_vc(vc, "Residual")
                total_var <- var_loc + var_resid
                if (total_var == 0) next
                results[i, "profile_location"] <- var_loc   / total_var
                results[i, "Residual"]         <- var_resid / total_var
            }
        }, error = function(e) {})
    }
    return(results)
}

# ── Main loop ─────────────────────────────────────────────────────────────────

out_dir <- "final_results/plot_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (cfg in data_configs) {
    cat("\n===", cfg$label, "===\n")

    ss_path <- file.path("data/FMC_processed", paste0(cfg$prefix, "_ss.csv"))
    ds_path <- file.path("data/FMC_processed", paste0(cfg$prefix, "_ds.csv"))

    if (!file.exists(ss_path) || !file.exists(ds_path)) {
        cat("Skipping", cfg$label, "(files not found)\n"); next
    }

    data <- rbind(read.csv(ss_path, na.strings = ""),
                  read.csv(ds_path, na.strings = ""))
    data <- data[!is.na(data$genus), ]
    data <- data[, setdiff(names(data), c("core_replicate", "tech_replicate", "library_type"))]
    data <- merge(data, rsa, by = "library_id", all.x = TRUE)
    data <- data[!is.na(data$profile_location), ]
    data <- data[data$profile_location != "54m", ]

    genus_totals <- tapply(data$n_reads, data$genus, sum, na.rm = TRUE)
    genus_n_libs <- tapply(data$library_id, data$genus, function(x) length(unique(x)))
    keep_genera  <- names(genus_totals)[genus_totals >= minReadsPerTaxon &
                                         genus_n_libs >= minSamples]
    data <- data[data$genus %in% keep_genera, ]

    count_wide <- data %>%
        group_by(library_id, genus) %>%
        summarise(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = library_id, values_from = n_reads,
                           values_fill = 0L)

    reads_table <- as.matrix(count_wide[, -1])
    rownames(reads_table) <- count_wide$genus
    reads_table <- reads_table[, colSums(reads_table) > 0, drop = FALSE]

    lib_meta <- data %>% distinct(library_id, library_type, profile_location, core_replicate)
    lib_meta  <- lib_meta[match(colnames(reads_table), lib_meta$library_id), ]

    cat("  Taxa:", nrow(reads_table), " | Libraries:", ncol(reads_table), "\n")
    print(table(lib_meta$library_type))

    var_df <- get_var_explained(reads_table, lib_meta)

    rds_path <- file.path(out_dir, paste0("var_partitioning_", cfg$label, ".rds"))
    saveRDS(list(var_df = var_df, label = cfg$label), rds_path)
    cat("  Saved:", rds_path, "\n")
}

cat("\nDone.\n")
