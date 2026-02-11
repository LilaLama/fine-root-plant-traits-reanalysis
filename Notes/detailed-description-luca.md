## Assessment of imputation uncertainty

In the original publication, the authors applied one single missForest imputation on data, that were missing not more than 2 below-ground and not more than 3 above-ground traits thereby increasing the number of species from 301 to 1,218. 
<how missForest works>
<why they chose missForest>

Even if imputation does not systematically distort ordination geometry, imputation uncertainty could still affect how stable individual species’ positions are in PCA space. To evaluate whether uncertainty in imputed values propagates into meaningful changes in PCA geometry, we explicitly tested ordination stability under different levels of imputation uncertainty.

The steps correspond to the R scripts:
- `code/061_multiple_imputation.R` — parametric bootstrap around missForest OOB error ("1× OOB")
- `code/062_imputation_mean_simple.R` — naive mean imputation baseline
- `code/063_multiple_imputation_2x_OOB.R` — sensitivity analysis with doubled OOB uncertainty ("2× OOB")

---

### Step 1 — Naive mean imputation baseline (Script 062)

**(i) Why do we do this step?**

- To check whether the global PCA structure (axes, trait loadings, species gradients) is robust to the choice of imputation method.
- Mean imputation is a deliberately naive baseline: it ignores phylogeny and multivariate structure and simply replaces each missing value by the trait-wise mean.
- If PCA results under mean imputation are very similar to those under missForest, this suggests that the ordination is mainly driven by the observed data structure and is not overly sensitive to sophisticated imputation.
- If results differ strongly, this indicates that missForest (and the information it uses, including phylogeny) meaningfully changes the functional space and therefore needs to be justified more carefully.

**(ii) How do we do it? (code and outputs)**

Script: `code/062_imputation_mean_simple.R`

1. Load the trait matrix used in the original PCA and the set of selected traits:

```r
traitsUse    <- readRDS("data/traitsUse.rds")
traitsSelect <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")
```

2. Compute trait-wise means from observed (non-missing) values and replace NAs with these means:

```r
trait_means <- sapply(traitsUse[, traitsSelect], function(col) mean(col, na.rm = TRUE))

imputedTraits_mean <- traitsUse[, traitsSelect]
for (j in colnames(imputedTraits_mean)) {
	imputedTraits_mean[is.na(imputedTraits_mean[, j]), j] <- trait_means[j]
}

AllTraitsAllInfo_mean <- AllTraitsAllInfo
AllTraitsAllInfo_mean[, traitsSelect] <- imputedTraits_mean
saveRDS(AllTraitsAllInfo_mean, file = "data/imputedTraits_mean.rds")
```

3. Re-run the full PCA pipeline used in the main analysis (number of components via `paran`, varimax rotation, score scaling, fixed orientation rules, 4D and 2D TPDs) on the mean-imputed dataset:

```r
AllTraits_mean <- AllTraitsAllInfo_mean[, traitsSelect]
PCATotal_mean <- list()
PCATotal_mean$traits     <- AllTraits_mean
PCATotal_mean$dimensions <- paran(PCATotal_mean$traits)$Retained
PCATotal_mean$PCA <- psych::principal(scale(PCATotal_mean$traits),
																			nfactors = PCATotal_mean$dimensions,
																			rotate   = "varimax", covar = TRUE)

# Orientation rules (same as original script)
for (i in 1:PCATotal_mean$dimensions) {
	if (i == 1 & PCATotal_mean$PCA$loadings["ph", i]  < 0) { ... }
	if (i == 2 & PCATotal_mean$PCA$loadings["sla", i] < 0) { ... }
	if (i == 3 & PCATotal_mean$PCA$loadings["SRL", i] > 0) { ... }
    ...
}

saveRDS(PCATotal_mean, file = "data/PCATotal_mean_imputation.rds")
```

**Key outputs**

- `data/imputedTraits_mean.rds` — mean-imputed trait matrix attached to species metadata.
- `data/PCATotal_mean_imputation.rds` — full PCA object (scores, loadings, TPDs) under mean imputation.

**(iii) What do we learn from this step?**

- The major axes of variation and trait loadings under mean imputation are very similar to those under missForest.
- Species positions in the first four PCA dimensions show only small shifts; global gradients and group separations are preserved.
- This indicates that the geometry of the functional trait space is robust to replacing missForest with a much simpler imputation, and that the main ordination-based conclusions do not depend critically on the specific imputation algorithm.

---

### Step 2 — Propagating OOB imputation uncertainty (1× OOB; Script 061)

**(i) Why do we do this step?**

- Single imputation (even with missForest) treats imputed values as if they were known, and therefore ignores imputation variance. This can make downstream uncertainty (e.g. confidence intervals, spread of species scores) appear too small.
- missForest internally reports out-of-bag (OOB) prediction error per trait (NRMSE), which quantifies how well the random forest can predict missing values.
- We use these OOB errors to parameterize a parametric bootstrap around the single imputed dataset: We generate many alternative versions of the dataset in which only originally missing entries are perturbed according to their estimated imputation error.
- By re-running the full PCA pipeline on each bootstrap dataset, we can quantify how much imputation uncertainty propagates into species scores and trait loadings.

**(ii) How do we do it? (code and outputs)**

Script: `code/061_multiple_imputation.R`

1. Run missForest once with `variablewise = TRUE` to obtain a single imputed dataset and OOB errors per trait:

```r
library(missForest)

traitsAux    <- readRDS("data/traitsAux.rds")
traitsSelect <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")

mf <- missForest(xmis = traitsAux, variablewise = TRUE, verbose = TRUE)
ximp0 <- mf$ximp
saveRDS(mf, file = "data/imputedTraitsOOB.rds")

oob_raw <- mf$OOBerror
names(oob_raw) <- colnames(ximp0)
oob_norm <- oob_raw[traitsSelect]  # extract OOB only for real traits
```

2. Convert trait-wise NRMSE to absolute imputation error SDs (on the log10 scale of the traits):

```r
sd_obs <- sapply(colnames(ximp0), function(j) sd(traitsAux[, j], na.rm = TRUE))
noise_sd <- oob_norm * sd_obs              # RMSE per trait
noise_sd[is.na(noise_sd) | noise_sd <= 0] <- 1e-8
missing_mat <- is.na(as.matrix(traitsAux)) # positions of originally missing values
```

3. Define a helper that re-runs the full PCA/TPD pipeline on an arbitrary imputed dataset (copied from the main script so that all steps are identical):

```r
make_PCAtotal_from_AllTraits <- function(AllTraitsAllInfo_obj, suffix = NULL) {
	# 1) extract traits, 2) run PCA (paran + principal with varimax),
	# 3) apply the same orientation rules as in the main analysis,
	# 4) compute 4D TPDs; return a PCATotal list and optionally save it
}
```

4. Bootstrap around the single imputation by adding parametric noise only to originally missing entries. For each of `nboot = 50` replicates:

```r
nboot  <- 50
outdir <- "data/imputed_bootstrap"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

scores_list   <- vector("list", nboot + 1)
loadings_list <- vector("list", nboot + 1)

# Single (deterministic) imputation as baseline
AllTraitsAllInfo_single <- AllTraitsAllInfo
AllTraitsAllInfo_single[, traitsSelect] <- ximp0[, traitsSelect]
saveRDS(AllTraitsAllInfo_single, file = file.path(outdir, "imputed_single.rds"))

PCA_single <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_single, suffix = "_single")
scores_list[[1]]   <- PCA_single$PCA$scores
loadings_list[[1]] <- PCA_single$PCA$loadings[, 1:4]

cat("\nBootstrapping imputations...\n")
pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
for (b in seq_len(nboot)) {
	utils::setTxtProgressBar(pb, b)
	set.seed(10000 + b)
	xboot <- ximp0
	for (j in seq_len(ncol(xboot))) {
		miss_idx <- which(missing_mat[, j])
		if (length(miss_idx) > 0) {
			xboot[miss_idx, j] <- rnorm(length(miss_idx),
										mean = ximp0[miss_idx, j],
										sd   = noise_sd[j])
		}
	}
	AllTraitsAllInfo_b <- AllTraitsAllInfo
	AllTraitsAllInfo_b[, traitsSelect] <- xboot[, traitsSelect]
	saveRDS(AllTraitsAllInfo_b,
					file = file.path(outdir, sprintf("imputed_boot_%03d.rds", b)))

	PCAb <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_b,
																			 suffix = sprintf("_boot_%03d", b))
	scores_list[[b + 1]]   <- PCAb$PCA$scores
	loadings_list[[b + 1]] <- PCAb$PCA$loadings[, 1:4]
}
close(pb)
```

5. Finally, we convert the lists of scores and loadings into arrays, compute the mean and standard deviation across all bootstrap replicates for each species × component and each trait × component, and save these summaries to `PCA_scores_boot_summary.rds` and `PCA_loadings_boot_summary.rds`.

**Key outputs**

- `data/imputedTraitsOOB.rds` — missForest object with single imputation and trait-wise OOB errors.
- `data/imputed_bootstrap/imputed_single.rds` — baseline single-imputed dataset.
- `data/imputed_bootstrap/imputed_boot_001.rds` … `imputed_boot_050.rds` — 50 bootstrap datasets with perturbed imputed values.
- `data/imputed_bootstrap/PCATotal_ImputedObs_single.rds` and `PCATotal_ImputedObs_boot_*.rds` — PCA results for single and bootstrap datasets.
- `data/imputed_bootstrap/PCA_scores_boot_summary.rds` — mean and SD of species scores (PC1–PC4) across bootstraps.
- `data/imputed_bootstrap/PCA_loadings_boot_summary.rds` / `.csv` — mean and SD of trait loadings across bootstraps.

**(iii) What do we learn from this step?**

- For all four main components, the bootstrap spread of species scores is small relative to the overall range of the ordination axes.
- Mean bootstrap scores are extremely close to the deterministic single-imputation scores, indicating that missForest does not systematically bias species positions.
- Trait loadings show very small bootstrap SDs; the identity and interpretation of the main PCA axes are highly stable under imputation uncertainty.
- Overall, this analysis shows that realistic levels of imputation uncertainty (as quantified by the OOB error) do not alter global gradients, trait axes, or group-level patterns.

---

### Step 3 — Sensitivity analysis with inflated uncertainty (2× OOB; Script 063)

**(i) Why do we do this step?**

- OOB error estimates are internal to missForest and could, in principle, underestimate true imputation uncertainty.
- To test the robustness of our conclusions to possible underestimation, we repeat the parametric bootstrap but inflate all OOB-based noise SDs by a factor of 2.
- This provides a conservative “worst-case” scenario: if PCA results remain stable even when imputation noise is doubled, then the main ordination conclusions are very robust.

**(ii) How do we do it? (code and outputs)**

Script: `code/063_multiple_imputation_2x_OOB.R`

1. As in Step 2, run missForest once, but scale the OOB errors by 2 before constructing the noise SDs:

```r
mf <- missForest(xmis = traitsAux, variablewise = TRUE, verbose = TRUE)
ximp0 <- mf$ximp

oob_raw <- mf$OOBerror * 2     # multiply OOB NRMSE by 2
names(oob_raw) <- colnames(ximp0)
oob_norm <- oob_raw[traitsSelect]

sd_obs <- sapply(colnames(ximp0), function(j) {
	v <- traitsAux[, j][!is.na(traitsAux[, j])]
	if (length(v) > 1) sd(v) else 0
})

noise_sd <- oob_norm * sd_obs  # 2× amplified RMSE per trait
noise_sd[is.na(noise_sd) | noise_sd <= 0] <- 1e-8
missing_mat <- is.na(as.matrix(traitsAux))

cat("\n*** OOB-BASED NOISE MULTIPLIED BY 2.0 ***\n")
print(noise_sd[traitsSelect])
```

2. Repeat the bootstrap loop and PCA pipeline as in Step 2

**Key outputs**

- `data/imputedTraits_2xOOB.rds` — missForest object used for the 2× OOB analysis.
- `data/imputed_bootstrap_2x/` — complete set of bootstrap imputed datasets and PCA results under 2× OOB noise.
- `data/imputed_bootstrap_2x/PCA_scores_boot_summary.rds` and `PCA_loadings_boot_summary.rds` — summaries of species scores and trait loadings under 2× OOB.
- `data/imputed_bootstrap_2x/PCA_scores_1x_vs_2x_boot.rds` — direct numeric comparison of 1× vs 2× OOB bootstrap means and SDs.
- `data/imputed_bootstrap_2x/PCA_scores_1x_vs_2x_*.png` — plots of 1× vs 2× OOB bootstrap mean scores with error bars.

**(iii) What do we learn from this step?**

- Doubling the OOB-based noise modestly increases the SD of species scores, but the increase is small relative to the total spread of the PCA axes.
- Mean species positions under 2× OOB remain extremely close to those under 1× OOB; the `delta` and `rel_ratio` summaries show only minor shifts.
- Trait loadings under 2× OOB remain highly similar to those under 1× OOB, and the qualitative interpretation of PCA axes is unchanged.
- The overall geometry of the ordination (orientation of trait axes, position of species groups, large-scale gradients) is therefore robust even to deliberately inflated imputation uncertainty.

---

### Overall conclusion from the imputation uncertainty assessment

Across the three steps (mean baseline, 1× OOB bootstrap, and 2× OOB sensitivity), the re-analysis shows that:

- The global PCA structure is insensitive to the choice between naive mean imputation and missForest.
- When imputation uncertainty is explicitly propagated using missForest OOB errors, species scores and trait loadings show only modest variation around the deterministic single-imputation solution.
- Even when OOB-based uncertainty is doubled, ordination geometry and trait axes remain highly stable.

Together, these results support the conclusion that the main findings about the fine-root functional trait space are robust to missing-data treatment and that imputation uncertainty primarily adds small, local noise rather than changing the underlying ecological signal.

