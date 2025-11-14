# Methods

## Data input
- Trait sources: TRY (aboveground traits) and GRoot (fine‑root traits).
- Taxonomy: standardized using The Plant List (Taxonstand) to harmonize species names.
- Occurrence/biome assignment: species assigned to biomes using GBIF occurrence records (see script for exact filtering/thresholds).

## Preprocessing
- Log10‑transform: applied to skewed continuous traits (improves normality and multiplicative-scale interpretation).
- Centering & scaling: traits centered (mean = 0) and scaled (sd = 1) before PCA so axes reflect standardized trait variation.
- Species filtering: analyses use subsets meeting minimum completeness thresholds (e.g., ≥50% above / below traits).
- Imputation: single missForest imputation used to fill missing trait cells. Phylogenetic information is included as PCoA axes (phylo eigenvectors) appended to predictors.

## Creation of Global Trait Space
- PCA: principal component analysis performed on standardized traits (psych::principal).
- Rotation: Varimax rotation applied to retained components to facilitate interpretation (redistributes variance among components).
- Comparison of spaces: seperate PCAs for subsets (above-only, below-only), compare species positions to full trait space using Procrustes tests (vegan::protest).
- TPD (Trait Probability Density): multivariate kernel density estimation. TPD = describes the probability that a randomly chosen species occurence falls at each point of the trait space. Used to quantify overlap between species and to show where they typically sit.
- Multivariate nulls: simulate many datasets from a simple multivariate normal model with the same overall mean and covariance as the real data; compare the observed trait patterns (e.g. PCA shape or variance) to those simulated to test if the real data differ from random multivariate structure.

## Exploration of patterns / downstream analyses
1. PERMANOVA (vegan::adonis) on Euclidean distances from trait coordinates to quantify how much variation in trait space is explained by growth form, family, biomes.
2. Dissimilarity among groups: compute community TPDs (sum of species' trait probability kernels) and pairwise TPD dissimilarities between groups; compare these trait dissimilarities to climate/abiotic dissimilarities (e.g., Gower distance) to test whether different environments host different trait pools.
3. Redundancy (TPDRichness): compute trait‑space volume occupied by each group (TPDRichness). Compare observed volumes to a null distribution from randomized species assignments to quantify how much trait space is unique (low redundancy) versus overlapped by many species (high redundancy).

## Uncertainty propagation (imputation)
- Single imputation treats imputed cells as known; to propagate uncertainty:
  1. Run missForest once, extract mf$OOBerror (per‑variable normalized RMSE).
  2. Derive absolute noise SD per trait: noise_sd = oob_norm * sd(observed trait).
  3. For nboot replicates, perturb only originally missing cells by sampling Normal(mean = imputed value, sd = noise_sd) → create multiple plausible imputed datasets.
  4. Run PCA and downstream analyses on each replicate; aggregate means and SDs of PCA scores, FRic, dissimilarities to quantify uncertainty.
- Practical note: set seed only if bit‑exact reproducibility is required; otherwise save mf (saveRDS) and the bootstrap outputs.