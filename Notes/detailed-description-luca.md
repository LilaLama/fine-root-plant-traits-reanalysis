# Assessment of imputation uncertainty

## Aim of imputation analysis

In the original publication, the authors applied one single missForest imputation on data, that were missing not more than 2 below-ground and not more than 3 above-ground traits thereby increasing the number of species from 301 to 1,218. 

The goal of the imputation analysis is to assess how sensitive our conclusions about the functional trait space are to the way missing data are handled and to the remaining uncertainty in imputed values.

More specifically, we aim to:
- test whether the geometry of the functional space (PCA axes, trait loadings, species gradients) is robust to replacing the original single missForest imputation with a deliberately naive mean-imputation baseline;
- quantify how missForest out-of-bag (OOB) imputation uncertainty propagates into species scores and trait loadings in PCA space; and
- evaluate a conservative worst-case scenario by doubling OOB-based uncertainty, to check whether our main ordination-based conclusions remain stable even under inflated imputation noise.

These aims are implemented in three complementary steps, corresponding to the following R scripts:
- `code/061_multiple_imputation.R` — parametric bootstrap around missForest OOB error ("1× OOB")
- `code/062_imputation_mean_simple.R` — naive mean imputation baseline
- `code/063_multiple_imputation_2x_OOB.R` — sensitivity analysis with doubled OOB uncertainty ("2× OOB")

---

## Out-of-bag (OOB) error in missForest

missForest (R package missForest; Stekhoven & Bühlmann 2012) is a random forest–based imputation algorithm. Like standard random forests, it uses bootstrap samples of the data to grow each tree: for every tree, a bootstrap sample of observations is drawn (with replacement) to build the tree, and the remaining observations are "out-of-bag" (OOB) for that tree. For an observation, an OOB prediction is obtained by aggregating predictions from all trees for which this observation was not used in tree construction. Comparing these OOB predictions to the observed values yields an internal cross-validation estimate of prediction error.

In missForest, this internal error is reported as an overall OOB error and, if `variablewise = TRUE`, as a separate error for each imputed variable. For continuous variables the error metric is the normalized root mean squared error (NRMSE), i.e. the root mean squared difference between OOB predictions and observed values, divided by the sample standard deviation of the observed values. For categorical variables missForest reports the proportion of falsely classified entries (PFC). In our case all imputed traits are continuous, so we use the trait-wise NRMSEs returned by missForest as measures of imputation uncertainty for each trait on the log10 scale.

In the OOB-bootstrap part of our analysis (Step 2 and 3), we convert these trait-wise NRMSEs back to absolute error scales by multiplying them with the observed trait standard deviation. This gives one noise standard deviation per trait (an estimated RMSE), which we then use to generate parametric noise around the single imputed values: for each originally missing entry of trait, we draw bootstrap values from a normal distribution centred on the missForest imputation with standard deviation equal to this trait-specific error. This procedure allows us to propagate the uncertainty quantified by the missForest OOB error through the PCA and all downstream summaries.

---

### Step 1 — Naive mean imputation baseline (Script 062)

**(i) Why do we do this step?**

- To check whether the global PCA structure (axes, trait loadings, species gradients) is robust to the choice of imputation method.
- Mean imputation is a deliberately naive baseline: it ignores phylogeny and multivariate structure and simply replaces each missing value by the trait-wise mean.
- If PCA results under mean imputation are very similar to those under missForest, this suggests that the ordination is mainly driven by the observed data structure and is not overly sensitive to sophisticated imputation.
- If results differ strongly, this indicates that missForest (and the information it uses, including phylogeny) meaningfully changes the functional space and therefore needs to be justified more carefully.

**(ii) How do we do it? (code and outputs)**

In `code/062_imputation_mean_simple.R`, we take the trait matrix and selected traits used in the original PCA, replace all missing entries by their trait-wise means computed from the observed values, and then re-run exactly the same PCA pipeline as in the main analysis (component selection with `paran`, varimax-rotated `psych::principal`, orientation rules, and saving the resulting PCA object). The core steps are sketched below:

```r
traitsUse        <- readRDS("data/traitsUse.rds")
traitsSelect     <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")

## Mean imputation
trait_means <- sapply(traitsUse[, traitsSelect], function(col) mean(col, na.rm = TRUE))
imputedTraits_mean <- traitsUse[, traitsSelect]
for (j in colnames(imputedTraits_mean)) {
  imputedTraits_mean[is.na(imputedTraits_mean[, j]), j] <- trait_means[j]
}
AllTraitsAllInfo_mean <- AllTraitsAllInfo
AllTraitsAllInfo_mean[, traitsSelect] <- imputedTraits_mean
saveRDS(AllTraitsAllInfo_mean, "data/imputedTraits_mean.rds")

## PCA on mean-imputed traits (same pipeline as main analysis)
AllTraits_mean            <- AllTraitsAllInfo_mean[, traitsSelect]
PCATotal_mean             <- list()
PCATotal_mean$traits      <- AllTraits_mean
PCATotal_mean$dimensions  <- paran(PCATotal_mean$traits)$Retained
PCATotal_mean$PCA <- psych::principal(scale(PCATotal_mean$traits),
                                      nfactors = PCATotal_mean$dimensions,
                                      rotate   = "varimax", covar = TRUE)

# Orientation rules and TPD calculations follow the main script
saveRDS(PCATotal_mean, "data/PCATotal_mean_imputation.rds")
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

- A single missForest imputation gives one "best guess" dataset and treats imputed values as fixed, so downstream analyses ignore imputation variance and may understate uncertainty.
- missForest provides trait-wise OOB error estimates that quantify the typical prediction error for each trait (see OOB section above). In this step we turn these error estimates into a parametric bootstrap around the single imputed dataset.
- By repeatedly perturbing only the originally missing entries according to their trait-specific OOB error and re-running the full PCA pipeline, we obtain a direct measure of how much imputation uncertainty propagates into species scores and trait loadings.

**(ii) How do we do it? (code and outputs)**

In `code/061_multiple_imputation.R`, we first run missForest with `variablewise = TRUE` to obtain a single imputed dataset together with trait-wise OOB NRMSE values. We then convert these normalized errors back to absolute error scales, identify which entries were originally missing, and, for each of 50 bootstrap replicates, add Gaussian noise with trait-specific SD only to those missing cells. For every perturbed dataset we re-run the full PCA/TPD pipeline using the same helper as in the main analysis, and finally summarize the resulting lists of scores and loadings across bootstraps (means and SDs) and save them for downstream use. The key implementation steps are:

```r
library(missForest)

## Single missForest imputation + OOB errors
traitsAux        <- readRDS("data/traitsAux.rds")
traitsSelect     <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")

mf    <- missForest(xmis = traitsAux, variablewise = TRUE, verbose = TRUE)
ximp0 <- mf$ximp
saveRDS(mf, "data/imputedTraitsOOB.rds")

oob_norm <- mf$OOBerror[traitsSelect]              # trait-wise NRMSE
sd_obs   <- sapply(traitsSelect, function(j) sd(traitsAux[, j], na.rm = TRUE))
noise_sd <- oob_norm * sd_obs                      # RMSE per trait
noise_sd[is.na(noise_sd) | noise_sd <= 0] <- 1e-8
missing_mat <- is.na(as.matrix(traitsAux))         # originally missing entries

## Parametric OOB bootstrap around the single imputation
nboot  <- 50
outdir <- "data/imputed_bootstrap"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

scores_list   <- vector("list", nboot + 1)
loadings_list <- vector("list", nboot + 1)

AllTraitsAllInfo_single <- AllTraitsAllInfo
AllTraitsAllInfo_single[, traitsSelect] <- ximp0[, traitsSelect]
saveRDS(AllTraitsAllInfo_single, file.path(outdir, "imputed_single.rds"))

PCA_single          <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_single, suffix = "_single")
scores_list[[1]]    <- PCA_single$PCA$scores
loadings_list[[1]]  <- PCA_single$PCA$loadings[, 1:4]

for (b in seq_len(nboot)) {
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
  saveRDS(AllTraitsAllInfo_b, file.path(outdir, sprintf("imputed_boot_%03d.rds", b)))

  PCAb                 <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_b, suffix = sprintf("_boot_%03d", b))
  scores_list[[b + 1]] <- PCAb$PCA$scores
  loadings_list[[b+1]] <- PCAb$PCA$loadings[, 1:4]
}

## Summaries across bootstraps (not shown in detail here)
saveRDS(list(scores = scores_list, loadings = loadings_list),
        file.path(outdir, "PCA_scores_loadings_boot_raw.rds"))
```

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

- The OOB-based error estimates used in Step 2 are model-based and might, in principle, underestimate the true uncertainty in imputed values.
- To probe this, we repeat exactly the same OOB-bootstrap procedure but inflate all trait-wise OOB error SDs by a factor of two before generating the perturbations.
- If PCA scores, loadings and group-level patterns remain essentially unchanged under this amplified noise, we gain evidence that our main ordination conclusions are robust even under a conservative "worst-case" imputation uncertainty.

**(ii) How do we do it? (code and outputs)**

In `code/063_multiple_imputation_2x_OOB.R`, we repeat the same workflow as in Step 2 but inflate the trait-wise OOB NRMSEs by a factor of two before converting them to absolute error SDs. The subsequent parametric bootstrap, PCA re-fitting and summarisation are identical to the 1× OOB case, so that the only difference between both analyses is the assumed size of the imputation noise. The essential code changes are:

```r
library(missForest)

traitsAux        <- readRDS("data/traitsAux.rds")
traitsSelect     <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")

mf    <- missForest(xmis = traitsAux, variablewise = TRUE, verbose = TRUE)
ximp0 <- mf$ximp

## Inflate trait-wise OOB errors by factor 2
oob_raw  <- mf$OOBerror * 2
oob_norm <- oob_raw[traitsSelect]

sd_obs   <- sapply(traitsSelect, function(j) {
  v <- traitsAux[, j][!is.na(traitsAux[, j])]
  if (length(v) > 1) sd(v) else 0
})

noise_sd <- oob_norm * sd_obs   # 2x amplified RMSE per trait
noise_sd[is.na(noise_sd) | noise_sd <= 0] <- 1e-8
missing_mat <- is.na(as.matrix(traitsAux))

cat("\n*** OOB-BASED NOISE MULTIPLIED BY 2.0 ***\n")
print(noise_sd[traitsSelect])

## The subsequent bootstrap loop, PCA re-fitting and summaries
## follow exactly the same structure as in code/061_multiple_imputation.R,
## but using these inflated noise_sd values.
```

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
- The exact ordinations obtained under 1× and 2× OOB will be compared formally in a later step using a Procrustes test.

---

### Overall conclusion from the imputation uncertainty assessment

Across the three steps (mean baseline, 1× OOB bootstrap, and 2× OOB sensitivity), the re-analysis shows that:

- The global PCA structure is insensitive to the choice between naive mean imputation and missForest.
- When imputation uncertainty is explicitly propagated using missForest OOB errors, species scores and trait loadings show only modest variation around the deterministic single-imputation solution.
- Even when OOB-based uncertainty is doubled, the visual ordination geometry and trait axes appear highly stable.

Taken together, these results strongly suggest that the main findings about the fine-root functional trait space are robust to missing-data treatment and that imputation uncertainty primarily adds small, local noise rather than changing the underlying ecological signal. A definitive quantitative assessment of how similar the ordinations are will, however, come from the Procrustes comparisons presented in the subsequent analysis.

