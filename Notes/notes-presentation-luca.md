## Imputation uncertainty block – talk notes

### Slide 1 – Assessing imputation uncertainty

Title: Assessing imputation uncertainty

Storyline:
- Original study: single missForest imputation → imputed values are treated like observed ones.
- That means: PCA scores, loadings and all ecological conclusions act as if there was no uncertainty in the missing traits.

Dataset overview (use figure):
- Complete dataset (no/mild missing data): 301 species (17.5%).
- Imputed dataset (≤50% missing traits, filled with missForest): 1,218 species (70.9%).
- Excluded species (>50% missing, not imputed): 200 species (11.6%).
- So most of the functional space is informed by imputed values, not by fully observed species.

Goal of this analysis step (say explicitly):
- Test how robust the functional trait space (PCA axes, trait loadings, species gradients) is to
	1) the **choice of imputation method**, and
	2) the **remaining uncertainty in imputed values**.

Transition: “To answer this, we designed a three-step sensitivity analysis.”

---

### Slide 2 – Analytical approach (3 steps)

Title: Three-step assessment of imputation uncertainty

Short bullets:
- Step 1 – “Naive” mean imputation
	- Replace all NAs with simple trait-wise means (no phylogeny, no multivariate structure), then re-run the full PCA pipeline.
	- Question: Does this naive baseline already give essentially the same functional space as missForest, or does missForest really matter?

- Step 2 – 1× OOB bootstrap (realistic uncertainty)
	- Run missForest once with `variablewise = TRUE`, use trait-wise OOB errors to add parametric noise only to originally missing entries, and re-fit the PCA for 50 bootstrap datasets.
	- Question: Given realistic imputation errors, how much do species positions and trait axes actually move in PCA space (mean vs. spread of scores/loadings)?

- Step 3 – 2× OOB “worst case”
	- Repeat the same bootstrap, but double all OOB-based noise SDs before perturbing the imputed values.
	- Question: Even under a conservative worst case (2× OOB), do the main ordination patterns and ecological conclusions still hold?

Transition: “Now I’ll show you what each of these steps does to the PCA.”

---

### Slide 3 – Results Step 1: mean imputation

Title: Step 1 – Naive mean imputation vs. missForest

What to show:
- Left/right (or overlaid) PCA of missForest vs. mean-imputed data (PC1–PC2, optionally PC1–PC3).
- Barplots or table of trait loadings PC1–PC4 for missForest vs. mean imputation.

Talking points:
- Visually, the ordinations are very similar: same main gradients, same separation of species groups; mean imputation shows slightly “blurrier” clustering, but no new axes appear.
- Trait loadings on PC1–PC4 differ only marginally between missForest and mean imputation; the ecological interpretation of each axis stays the same.
- Take-home: The geometry of the functional trait space is **robust to the choice of imputation method** – a very simple imputation already recovers the main structure.

Transition: “So the space is not sensitive to how exactly we fill in NAs. Next, we asked how much uncertainty remains *within* missForest, given its OOB error estimates.”

---

### Slide 4 – Results Step 2: 1× OOB bootstrap

Title: Step 2 – 1× OOB bootstrap (realistic uncertainty)

What to show:
- Scatterplot of deterministic missForest PCA scores (single imputation) with overlaid ‘clouds’ / ellipses or error bars from the 1× OOB bootstrap (PC1–PC2, PC3–PC4).
- Possibly a summary of SD of scores per PC (boxplot or table) and SD of trait loadings.

Talking points:
- The bootstrap clouds are very tight: for all four components, the SD of species scores is small relative to the total axis range.
- Mean bootstrap positions are almost identical to the single-imputation scores → missForest does not systematically shift species, it just adds small symmetric jitter.
- Trait loadings show very small bootstrap SDs; the identity and interpretation of PCA axes are extremely stable.
- Take-home: **Realistic imputation uncertainty has only minor effects** on the ordination – local noise, but no structural changes.

Transition: “But OOB errors might themselves be a bit optimistic. So we repeated the same analysis, but deliberately doubled the noise.”

---

### Slide 5 – Results Step 3: 2× OOB “worst case”

Title: Step 3 – 2× OOB bootstrap (conservative worst case)

What to show:
- Same type of plot as for Step 2, but now with 2× OOB clouds.
- Maybe a direct comparison of 1× vs. 2× SDs of scores (e.g. ratio or difference).

Talking points:
- As expected, the bootstrap SD of species scores increases, but even at 2× OOB it is still small compared to the total spread of the PCA axes.
- Mean species positions under 2× OOB remain extremely close to those under 1× OOB; groupings and main gradients are visually unchanged.
- Trait loadings remain highly similar; the qualitative interpretation of the axes does not change.
- Take-home: **Even under exaggerated imputation uncertainty, the ordination geometry is very stable.**

Transition to conclusion: “So across all three steps, the space is strikingly robust. The next question is: can we quantify ‘how similar’ these ordinations really are?”

---

### Slide 6 – Overall conclusion + handover

Title: Summary and quantitative comparison

Key messages to state clearly:
- The functional trait space is robust to (i) replacing missForest with naive mean imputation, and (ii) propagating realistic and even doubled imputation uncertainty.
- Imputation uncertainty mainly adds small, local noise to species positions, but does not create or destroy the main ecological gradients.

Lead-in to Procrustes analysis (handover to Charlie):
- “So far I have shown this robustness visually. To **quantify** how similar these ordinations really are, we compared them formally using Procrustes tests.”
- “Charlie will now explain how these Procrustes comparisons work and what they tell us about the similarity between the different imputation-based ordinations.”

