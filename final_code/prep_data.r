# prep_data.r
# Produce the four processed long-format CSVs from raw inputs only.
#
# Inputs (all read from their canonical paths):
#   data/FLB_combined_euk_highconf_plants.tsv
#   data/FLB_combined_mic_highconf_mics.tsv
#   data/FMC_processed/replicate_sample_assignments_newNames.csv
#
# Outputs (written to data/FMC_processed/):
#   ellesmere_long_ss.csv
#   ellesmere_long_ds.csv
#   ellesmere_mic_long_ss.csv
#   ellesmere_mic_long_ds.csv
#
# Column format: library_id, genus, core_replicate, tech_replicate, n_reads, profile_location
#   tech_replicate   = paste(core_replicate, tech_replicate_from_RSA, sep = "_")
#   profile_location is joined from the RSA at prep time and carried in the output CSV.
#
# Run from the ellesmere/ project root.

library(dplyr)
library(tidyr)

# ── Load RSA ──────────────────────────────────────────────────────────────────

rsa <- read.csv("data/FMC_processed/replicate_sample_assignments_newNames.csv",
                check.names = FALSE)
rsa <- rsa[, c("library_id", "library_type", "profile_location", "core_replicate", "tech_replicate")]

# ── Sanity checks on RSA ──────────────────────────────────────────────────────

dup_ids <- unique(rsa$library_id[duplicated(rsa$library_id)])
if (length(dup_ids) > 0)
    cat("WARNING: duplicate library_ids in RSA (keeping first):", paste(dup_ids, collapse = ", "), "\n")

for (col in c("library_id", "library_type", "profile_location", "core_replicate", "tech_replicate")) {
    n_na <- sum(is.na(rsa[[col]]))
    if (n_na > 0) cat("WARNING: RSA has", n_na, "NA(s) in column '", col, "'\n")
}

# library_id should start with FLB{profile_location}{core_replicate}_
expected_prefix <- paste0("FLB", rsa$profile_location, rsa$core_replicate, "_")
bad_prefix <- !startsWith(rsa$library_id, expected_prefix)
if (any(bad_prefix)) {
    cat("WARNING:", sum(bad_prefix), "library_id(s) don't start with FLB{location}{core_rep}_:\n")
    print(rsa[bad_prefix, c("library_id", "profile_location", "core_replicate")])
}

# library_type should appear immediately after the first underscore in library_id
after_underscore <- sub("^[^_]+_", "", rsa$library_id)
bad_lt <- !startsWith(after_underscore, rsa$library_type)
if (any(bad_lt)) {
    cat("WARNING:", sum(bad_lt), "library_id(s) don't encode library_type after '_':\n")
    print(rsa[bad_lt, c("library_id", "library_type")])
}

rsa <- rsa[!duplicated(rsa$library_id), ]
rsa$tech_rep_out <- paste(rsa$core_replicate, rsa$tech_replicate, sep = "_")

# ── Helper: TSV → 4 long CSVs ────────────────────────────────────────────────

process_tsv <- function(tsv_path, suffix, out_prefix) {
    cat("Reading", tsv_path, "\n")
    raw <- read.delim(tsv_path, check.names = FALSE, stringsAsFactors = FALSE)

    # Keep only TotalReads columns (one per library) plus TaxName
    reads_cols <- grep(paste0("_", suffix, "_TotalReads$"), names(raw), value = TRUE)
    if (length(reads_cols) == 0)
        stop("No TotalReads columns found matching suffix '", suffix, "' in ", tsv_path)

    raw <- raw[, c("TaxName", reads_cols)]

    # Parse library_id from column names
    pattern <- paste0("_", suffix, "_TotalReads$")
    lib_ids  <- sub(pattern, "", reads_cols)

    # Keep only libraries in RSA
    keep <- lib_ids %in% rsa$library_id
    if (!any(keep)) stop("No library_ids from TSV match RSA for suffix '", suffix, "'")
    not_in_rsa <- lib_ids[!keep]
    if (length(not_in_rsa) > 0)
        cat("  NOTE:", length(not_in_rsa), "TSV library_id(s) not in RSA (dropped):",
            paste(not_in_rsa, collapse = ", "), "\n")
    raw <- raw[, c("TaxName", reads_cols[keep]), drop = FALSE]
    lib_ids <- lib_ids[keep]

    # Rename columns to library_id for clean pivot
    names(raw)[names(raw) != "TaxName"] <- lib_ids

    # Pivot to long
    long <- pivot_longer(raw, cols = all_of(lib_ids),
                         names_to = "library_id", values_to = "n_reads")
    names(long)[names(long) == "TaxName"] <- "genus"

    # Join RSA metadata
    long <- left_join(long, rsa, by = "library_id")
    long <- long[!is.na(long$library_type), ]

    # Check for unexpected NAs in profile_location after join
    na_loc <- unique(long$library_id[is.na(long$profile_location)])
    if (length(na_loc) > 0)
        cat("  WARNING: NA profile_location after join for library_id(s):",
            paste(na_loc, collapse = ", "), "\n")

    # Build output columns
    long$tech_replicate <- long$tech_rep_out

    out_cols <- c("library_id", "genus", "core_replicate", "tech_replicate", "n_reads", "profile_location")

    # Split by library_type and write
    for (lt in c("ss", "ds")) {
        sub <- long[long$library_type == lt, out_cols]
        out_path <- file.path("data/FMC_processed",
                              paste0(out_prefix, "_", lt, ".csv"))
        write.csv(sub, out_path, row.names = FALSE)
        cat("  Wrote", nrow(sub), "rows →", out_path, "\n")
    }
}

# ── Run for both organism types ────────────────────────────────────────────────

process_tsv("data/FLB_combined_euk_highconf_plants.tsv", "euk", "ellesmere_long")
process_tsv("data/FLB_combined_mic_highconf_mics.tsv",   "mic", "ellesmere_mic_long")

cat("\nDone.\n")
