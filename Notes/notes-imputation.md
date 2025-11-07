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
   - Imputation reliably reproduces species’ positions even if up to 50% of traits are missing.  
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

Practical, minimal options:
- Multiple imputation (recommended): run missForest m times (e.g., m = 10–20), compute PCA on each imputed dataset, align PC axes across imputations, then summarize species scores and all metrics (mean ± SD or 95% CI). Combine scalar estimates with Rubin’s rules where applicable.
- Bootstrap + imputation: bootstrap rows, impute each sample, compute PCA/metrics, summarize distribution.
- Model-based uncertainty: use Bayesian/EM PCA variants that return uncertainty for scores/loadings (BPCA, pcaMethods, or Bayesian factor models).

How to do multiple-imputation + PCA in practice (R sketch):
```r
# minimal sketch (no filepath)
library(missForest); library(psych); library(vegan)

m <- 20
scores_list <- vector("list", m)

# reference: do one imputation to set orientation
imp0 <- missForest(traitsAux)$ximp[, traitsSelect]
pc0 <- psych::principal(scale(imp0), nfactors = 4, rotate = "varimax")
ref_scores <- pc0$scores

for(i in 1:m){
  imp <- missForest(traitsAux)$ximp[, traitsSelect]
  pc <- psych::principal(scale(imp), nfactors = 4, rotate = "varimax")
  scores_list[[i]] <- pc$scores
}

# align each imputed PCA to reference via Procrustes
aligned <- lapply(scores_list, function(s){
  pr <- vegan::procrustes(ref_scores, s)
  pr$Yrot
})

# build array: species x axes x m, then summarize
scores_array <- array(unlist(aligned), c(nrow(ref_scores), ncol(ref_scores), m))
mean_scores <- apply(scores_array, c(1,2), mean)
sd_scores   <- apply(scores_array, c(1,2), sd)

# downstream: for each imputation compute TPDs/FRic/dissimilarities,
# then report mean ± CI across m runs
```