# Final Analysis: Code Files and Plots

All scripts run from the `ellesmere/` project root.
Sample metadata always loaded from `data/FMC_processed/replicate_sample_assignments_newNames.csv`.
Eukaryote and microbe analyses are separate throughout. 
Analyses use both ss and ds libraries as replicates.

---

## Code Files

### `prep_data.r`
**What it does:** Reads the two raw TSV files and the sample metadata CSV, then produces
the four processed long-format CSVs used by all downstream analysis scripts. Each output
row is one taxon × one library, with columns: library_id, genus, core_replicate,
tech_replicate, n_reads, profile_location.

**Inputs:**
- `data/FLB_combined_euk_highconf_plants.tsv`
- `data/FLB_combined_mic_highconf_mics.tsv`
- `data/FMC_processed/replicate_sample_assignments_newNames.csv`

**Outputs:**
- `data/FMC_processed/ellesmere_long_ss.csv`
- `data/FMC_processed/ellesmere_long_ds.csv`
- `data/FMC_processed/ellesmere_mic_long_ss.csv`
- `data/FMC_processed/ellesmere_mic_long_ds.csv`

**Must run before:** all other scripts.

---

### `run_variance_partitioning.r`
**What it does:** Fits a per-taxon linear mixed model with ss and ds combined:
`read_proportions ~ library_type + (1|profile_location) + (1|profile_location:core_replicate)`.
The library_type component (ds vs ss) is a fixed effect; variance explained by it
is computed as `var(X*beta)` across observations. profile_location captures between-location
variance (depth/position in the sediment profile). profile_location:core_replicate captures
variance between core replicates within a location. Residual captures within-core-rep
(tech-rep / lane) variance.
Filters: ≥ 1000 reads per taxon across all libraries, present in ≥ 5 libraries.
Normalization: proportionsCPM.

**Runs for:** Euk (`ellesmere_long`) and Mic (`ellesmere_mic_long`).

**Outputs:**
- `final_results/plot_data/var_partitioning_Euk.rds`
- `final_results/plot_data/var_partitioning_Mic.rds`

**Must run before:** `plot_varpart_violin.r`, `plot_varpart_barplot.r`

---

### `plot_composition_errorbars.r`
**What it does:** Makes stacked bar composition plots (% reads by taxon) across profile
heights, with ss and ds libraries pooled. Reads are summed across both library types within
each (profile_location, core_replicate), so each core replicate contributes one combined
profile. Proportions are computed per pooled replicate, then mean and 95% CI of the mean
are taken across core replicates at each profile height. Error bars shown for the 10 most
abundant taxa.

**Runs for:** Euk and Mic.

**Outputs:**
- `final_plots/composition_errorbars_euk.pdf`
- `final_plots/composition_errorbars_mic.pdf`

---

### `run_aitchison_distances.r`
**What it does:** Computes Aitchison distances (Euclidean on CLR-transformed proportions,
+0.5 pseudocount) between all pairs of libraries, classifies pairs as tech-rep (same
profile_location + same core_replicate), core-rep (same profile_location + different
core_replicate), or non-rep (different profile_location). Cross-type (ss vs ds) pairs
within the same profile_location are excluded. Runs permutation tests (9999 permutations)
for all pairwise group comparisons. Saves all plot-ready data (summary stats and
significance bracket positions) to RDS so the plot can be re-made without repeating
the permutation tests.

**Runs for:** Euk and Mic.

**Outputs:**
- `final_results/plot_data/aitchison_distances.rds` — plot-ready data (summary + sig brackets)
- `final_results/aitchison_summary.csv` — mean distance, SE, n per group

**Must run before:** `plot_aitchison_distances.r`

---

### `plot_aitchison_distances.r`
**What it does:** Reads pre-computed data from `run_aitchison_distances.r` and renders
the grouped bar plot of mean Aitchison distance ± SE with significance brackets.
Faceted by ss+ds combined / ss only / ds only. Can be re-run to adjust aesthetics
without repeating the permutation tests.

**Requires:** `final_results/plot_data/aitchison_distances.rds`

**Output:**
- `final_plots/aitchison_distances.pdf`

---

### `plot_varpart_violin.r`
**What it does:** Reads the variance partitioning RDS files from `run_variance_partitioning.r`
and makes violin + box plots of the per-taxon variance components:
ds vs ss / Profile location / Core rep / Residual. Euk and Mic shown as side-by-side facets.

**Requires:** RDS files from `run_variance_partitioning.r`

**Output:**
- `final_plots/varpart_violin.pdf`

---

### `plot_varpart_barplot.r`
**What it does:** Reads the variance partitioning RDS files and makes a horizontal
stacked bar chart with one bar per taxon, showing the proportional contribution of
each variance component (ds vs ss / Profile location / Core rep / Residual). Taxa sorted by
Profile location variance (ascending), so taxa with the most spatially structured signal
appear at the top. One plot per organism type.

**Requires:** RDS files from `run_variance_partitioning.r`

**Outputs:**
- `final_plots/varpart_barplot_euk.pdf`
- `final_plots/varpart_barplot_mic.pdf`

---

### `run_ss_vs_ds.r`
**What it does:** Tests whether ss and ds libraries differ in composition after
controlling for profile_location. Reads are aggregated to
(profile_location, library_type, core_replicate) level. Aitchison distances are
computed between three pair types:
- `ss-ss`: ssA vs ssB — core-rep variation within ss
- `ds-ds`: dsA vs dsB — core-rep variation within ds
- `matched-cross`: ssX vs dsX — library-type distance for the same core-rep

Only locations with ≥ 2 core_replicate labels present in both ss and ds are included.
Statistical test: `lmer(distance ~ pair_type + (1|profile_location))`, likelihood ratio ANOVA.
Pairwise contrasts computed from REML model. If matched-cross ≈ ss-ss and ds-ds,
ss and ds are not systematically more different from each other than core replicates
within the same type.

**Runs for:** Euk and Mic.

**Outputs:**
- `final_results/ss_vs_ds_contrasts.csv` — pairwise contrast table (estimate, SE, t, p)
- `final_results/ss_vs_ds_means.csv` — mean Aitchison distance per pair type
- `final_results/plot_data/ss_vs_ds.rds` — raw pairs data frame + model summaries for re-plotting

---

## Plots

### `final_plots/composition_errorbars_euk.pdf`
Stacked bar chart showing eukaryote community composition (% reads) across profile heights,
with ss and ds libraries pooled. Each bar is a profile height; each color is a taxon
(plant/macroalgal genus). Error bars (95% CI of the mean) shown for the 10 most
abundant taxa. Illustrates how community composition changes along the sediment profile and
how consistent that composition is across core replicates.

### `final_plots/composition_errorbars_mic.pdf`
Same as above but for microbial families. Shows profile-height-structured turnover in
microbial community composition, with replication precision indicated by error bar width.

### `final_plots/aitchison_distances.pdf`
Grouped bar plot of mean Aitchison distance ± SE between tech-rep, core-rep, and non-rep
pairs, for Euk and Mic separately, faceted by ss+ds combined / ss only / ds only.
Significance brackets from 9999-permutation tests. Demonstrates that replicates (especially
tech-reps) are substantially more similar to each other than to samples from different
profile locations, validating the replication strategy. Also shows whether ss and ds
libraries behave similarly.

### `final_plots/varpart_violin.pdf`
Violin + box plots showing the distribution of per-taxon variance explained by each
component of the combined model: ds vs ss (library type), Profile location (spatial
position in the sediment profile), Core rep (A vs B within a location), and Residual.
Euk and Mic shown as side-by-side facets. Illustrates the relative magnitudes of different
sources of variation and how small the library-type effect is relative to spatial structure.

### `final_plots/varpart_barplot_euk.pdf`
Horizontal stacked bar chart, one bar per eukaryote taxon, showing each taxon's
variance partitioned into ds vs ss / Profile location / Core rep / Residual.
Sorted by Profile location variance. Shows which specific taxa are more vs less spatially
structured and which carry a detectable library-type signal.

### `final_plots/varpart_barplot_mic.pdf`
Same as above but for microbial families.
