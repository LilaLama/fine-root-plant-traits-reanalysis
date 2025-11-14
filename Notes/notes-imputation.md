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

missForest gives a best guess but also reports OOB error estimates for each imputed value. This uncertainty can be propagated through PCA and downstream analyses using multiple imputation (e.g., Rubin’s rules).

## How to propagate OOB uncertainty to PCA?
set seed? 
1. run missForest once
2. for each imputed value, sample from a distribution centered on the imputed value with spread given by the OOB error (e.g., normal distribution with mean = imputed value, sd = oob_norm * sd(imputed values for that trait))
3. on the missing values, for nboot replicates, generate nboot datasets by sampling imputed values as above
4. for each bootstrapped dataset run PCA (and downstream analyses)
5. get mean and sd of PCA scores over bootstraps to get uncertainty in species positions

