# Imputation — Justification for Single missForest Imputation (from peer review response)
1. **Purpose of the imputation:**  
   - The goal was **not to accurately reconstruct every single trait value**, but to determine each species’ **position in the functional space** (PCA).  
   - Even if some traits are missing, the combination of traits and phylogenetic information allows accurate placement of species in the space.  

2. **Simulation to test accuracy:**  
   - They artificially introduced missing values in species with complete trait data, mimicking the pattern in the original dataset.  
   - The missForest imputation was applied, and the species were projected into the same functional space as in the original dataset (no new PCA).  
   - They compared positions before and after imputation using **normalized root mean square error (NRMSE)**.  

3. **Results:**  
   - Very low errors:  
     - Components 1–3: 0.5–0.6% of the dimension range  
     - Component 4: 3.4%  
   - Imputation reliably reproduces species’ positions even if up to 50% of trait values are missing.  
   - Comparison to literature: Penone et al. (2014) reported NRMSE >6%, so this is much lower.  

4. **Phylogenetic information:**  
   - Imputation performs better for families with many species, but even single-species families show very low errors.  

**Conclusion:**  
- A single missForest imputation is sufficient because the goal is the **species’ positions in PCA space**, not exact trait reconstruction.  
- Simulations confirm that the results are robust, justifying a single imputation run.

# BUT
Ideally you should propagate imputation uncertainty into the PCA and all downstream metrics. Single missForest imputation ignores imputation variance and can make confidence intervals/p‑values too narrow, even if the simulated NRMSE is small.

Why it matters (very short)

Single imputation gives one “best guess” dataset → downstream results treat imputed values as known.
Proper propagation shows how much uncertainty in species positions / FRic / dissimilarities comes from missing data.

missForest gives a best guess but also reports OOB error estimates for each imputed value. This uncertainty can be propagated through PCA and downstream analyses using multiple imputation (e.g., Rubin’s rules).

## How to propagate OOB uncertainty to PCA?
set seed? 
1. run missForest once
2. for each imputed value, sample from a distribution centered on the imputed value with spread given by the OOB error (e.g., normal distribution with mean = imputed value, sd = oob_norm * sd(imputed values for that trait))
3. on the missing values, for nboot replicates, generate nboot datasets by sampling imputed values as above
4. for each bootstrapped dataset run PCA (and downstream analyses)
5. get mean and sd of PCA scores over bootstraps to get uncertainty in species positions

## Implementation: Three-method comparison (Scripts 061–063, 10a)

To address the uncertainty propagation question and test whether the OOB estimates are reliable, three imputation strategies were implemented:

### 1. **Original Single missForest (Script 06)**
   - Single deterministic imputation
   - PCA run once on this imputed dataset
   - Baseline for comparison

### 2. **Mean Imputation (Script 062)**
   - Naive baseline: replace all NAs with trait-wise means
   - Full PCA pipeline (varimax rotation, TPDs, etc.)
   - **Purpose:** Test whether missForest imputation is *essential* for PCA structure, or if naive mean-filling gives similar results
   - If loadings differ significantly from missForest → missForest provides meaningful information
   - If loadings are similar → PCA structure is robust to imputation method (reassuring)

### 3. **OOB Bootstrap Uncertainty Propagation (Scripts 061 & 063)**

#### Script 061 (1x OOB):
   - Run missForest once → get imputed values + OOB error estimates per variable
   - For 50 bootstrap replicates:
     - For each originally missing value, sample from $N(\text{imputed value}, \text{OOB error} \times \text{sd}_{\text{trait}})$
     - Run full PCA on each bootstrap dataset
   - Compute mean and SD of PCA scores across 50 bootstraps
   - Visualizations: Deterministic vs. bootstrap mean positions with ±1 SD error bars
   - Files saved: `PCA_scores_single_vs_boot.rds`, `PCATotal_ImputedObs_boot_*.rds`

#### Script 063 (2x OOB - Sensitivity analysis):
   - **Same as 063 but with OOB noise multiplied by 2.0**
   - Purpose: Test whether OOB uncertainty estimates are underestimated
   - If 2x OOB gives similar PCA results → OOB estimates are reasonable
   - If 2x OOB gives substantially different results → OOB might be underestimated
   - Files saved: `PCA_scores_1x_vs_2x_boot.rds`, comparison plots in `data/imputed_bootstrap_2x/`

### 4. **Loadings Comparison (Script 10a)**
   - Extracts varimax-rotated PCA loadings for all 10 traits across PC1–PC4
   - Compares loadings across all four methods:
     - Original missForest
     - Mean imputation
     - 1x OOB bootstrap (mean)
     - 2x OOB bootstrap (mean)
   - Output: `data/Loadings_all_PCs_comparison.csv`
   - Summary statistics: Mean absolute differences between methods
   - Shows whether imputation method affects **PCA structure** (loadings) or only **species positions** (scores)

## Interpretation for Peer Review

### If mean imputation ≈ missForest:
- missForest adds little value → simple mean imputation is sufficient
- Suggests PCA structure is **robust** to imputation approach

### If mean imputation ≠ missForest:
- missForest provides meaningfully different imputations → phylogenetic info helps
- Justifies using missForest over naive approaches

### If 1x OOB ≈ 2x OOB:
- OOB uncertainty estimates are **not underestimated**
- Bootstrap procedure is adequate for propagating uncertainty
- Gives confidence in reported PCA score uncertainties

### If 1x OOB ≠ 2x OOB:
- OOB might be underestimated
- Recommend reporting wider confidence intervals
- Or use more conservative 2x OOB as sensitivity check

