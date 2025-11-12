# Plant Trait Analysis Pipeline Overview

## Overall Pipeline Summary

This project analyzes the functional diversity of vascular plant traits across above- and belowground compartments. The pipeline consists of 10 main R scripts (00-09) plus auxiliary functions, investigating whether trait diversity is globally higher above- than belowground through functional space analysis.

**Data**: Plant trait data from TRY database (aboveground) and GROOT database (fine-root traits)

**Main Approach**: 
- Principal Component Analysis (PCA) with varimax rotation to create functional spaces
- Trait Probability Density (TPD) estimation for functional diversity analysis
- Comparison of aboveground vs belowground trait spaces
- Statistical testing with null models and multivariate analyses

---

## Script Descriptions

### 00_All analyses.R
- **Master script** that executes all analyses in sequence
- Sets up working directory and loads required packages
- Provides execution timing estimates for each step
- Creates output folder structure

### 01_Full empirical information spectrum.R
- **Data preparation and functional space creation** for species with complete trait information
- Loads trait data: aboveground traits (TRY), fine-root traits (GROOT), taxonomy
- Filters species with complete information: aboveground (2,630 spp), belowground (748 spp), both (301 spp)
- Performs PCA with varimax rotation for each trait subset
- Estimates TPD functions for functional diversity calculations
- Creates 2D projections for visualization

### 01b_Eigenanalysis from correlation matrix.R
- **Alternative PCA approach** using correlation matrix eigenanalysis
- Performs spectral decomposition on trait correlation matrix
- Applies varimax rotation to first 4 components
- Provides comparison to standard PCA results

### 02_Spectra Only above-Only below.R
- **Visualization of separate functional spaces** (aboveground-only, belowground-only)
- Creates trait probability density plots for each compartment
- Generates contour plots showing species density distributions
- Produces figures analogous to Díaz et al. (2016) and Bergmann et al. (2020)

### 03_4D Spectrum.R
- **Visualization of combined 4D functional space** (301 species with complete data)
- Plots all pairwise combinations of the 4 principal components
- Creates species density contour plots for each 2D plane
- Generates Extended Data figures and main manuscript Figure 1

### 04_Procrustes_Analyses.R
- **Comparison between functional spaces** using Procrustes analysis
- Tests similarity between aboveground-only vs full spectrum (C1-C2 plane)
- Tests similarity between belowground-only vs full spectrum (C3-C4 plane)
- Compares aboveground-only vs belowground-only spaces
- Quantifies congruence between different trait space approaches

### 05_Comparing_with_multivariate_Normal.R
- **Null model testing** against multivariate normal distributions
- Generates 499 random communities with same means/covariances as observed data
- Compares functional evenness (FEve) and functional divergence (FDiv)
- Tests whether observed trait distributions deviate from normal expectations
- **Warning**: Very computationally intensive (22.8+ hours)

### 06_Imputed information spectrum.R
- **Trait imputation and expanded analysis** (1,218 species)
- Selects species with ≥50% trait completeness for each compartment
- Uses phylogenetic information and missForest for trait imputation
- Creates functional space based on imputed dataset
- Provides larger sample size for subsequent analyses

### 07_Permanova_Analyses.R
- **PERMANOVA analysis** of trait variation explained by ecological factors
- Tests variance explained by: growth form (woody/herbaceous), plant family, biome
- Analyzes individual traits, principal components, and functional planes
- Uses 500 randomizations for biome assignments (species can occur in multiple biomes)

### 08_Dissimilarity_Analyses.R
- **Functional dissimilarity analysis** between ecological groups
- Calculates dissimilarity between woody vs herbaceous species
- Analyzes dissimilarity patterns among plant families (≥15 species)
- Examines biome functional dissimilarity and climate relationships
- Performs Mantel tests and creates Figure 2
- Tests correlation between trait dissimilarity and climate dissimilarity

### 09_Redundancy_Analyses.R
- **Functional redundancy analysis** using null models
- Estimates functional richness (FRic) for different ecological groups
- Compares observed vs expected functional space occupation
- Tests 4,999 random assemblages for statistical significance
- Analyzes redundancy patterns for growth forms, families, and biomes
- Generates Figure 3 and tests for functional clustering vs overdispersion

### Aux_Functions.R
- **Supporting functions** for TPD analysis and visualization
- Custom functions for handling large datasets (TPD_large, TPDc_large, etc.)
- Image plotting functions for functional space visualization
- Functional richness calculation functions
- Density profile functions for quantile analysis
- Modified TPD functions optimized for this specific analysis

---

## Key Outputs

- **Functional spaces**: PCA-based trait spaces for visualization and analysis
- **Statistical tests**: Procrustes, PERMANOVA, Mantel tests
- **Null models**: Comparison with random expectations
- **Figures**: Main manuscript figures and extended data visualizations
- **Data files**: Processed datasets saved as .rds files for reuse

---

## Critical Issue: Imputation Uncertainty Not Propagated

### Current Approach (Script 06)

The study uses a single `missForest` imputation run to handle missing trait data:

```r
# Current implementation - SINGLE imputation
phylDissAux <- sqrt(cophenetic(phylogenyAux))
pcoaPhyl <- cmdscale(phylDissAux, k=10) 
traitsAux <- cbind(traitsUse, pcoaPhyl)
colnames(traitsAux) <- c(colnames(traitsUse), paste0("PC.", 1:ncol(pcoaPhyl)))
traitsAux <- traitsAux[rownames(AllTraitsAllInfo), ]
imputedTraits <- (missForest(xmis= traitsAux)$ximp[, traitsSelect])
identical(rownames(AllTraitsAllInfo), rownames(imputedTraits))
AllTraitsAllInfo[, traitsSelect] <- imputedTraits # Single point estimates used
```

### Problems with Current Approach

1. **No uncertainty quantification**: Single imputation treats imputed values as known truth
2. **Overstated precision**: Downstream analyses ignore imputation variance
3. **Biased inference**: Standard errors and confidence intervals are too narrow
4. **Missing sensitivity analysis**: No assessment of how imputation affects conclusions

### Recommended Improvements

#### Option 1: Multiple Imputation with MICE-style Approach

```r
library(mice)
library(missForest)

# Multiple imputation approach
n_imputations <- 20  # Multiple datasets
imputation_results <- list()

# Create multiple imputed datasets
for(i in 1:n_imputations) {
  cat("Imputation", i, "of", n_imputations, "\n")
  
  # Use different random seeds for variability
  set.seed(1000 + i)
  
  # MissForest with different parameters
  imputed_data <- missForest(
    xmis = traitsAux,
    maxiter = 10,
    ntree = 100,
    variablewise = FALSE,
    verbose = FALSE
  )$ximp[, traitsSelect]
  
  imputation_results[[i]] <- imputed_data
}

# Pool results across imputations
pooled_results <- pool_functional_analyses(imputation_results)
```

#### Option 2: Bootstrap Imputation Uncertainty

```r
# Bootstrap-based uncertainty estimation
bootstrap_imputation <- function(data, n_bootstrap = 100) {
  
  bootstrap_results <- list()
  n_species <- nrow(data)
  
  for(b in 1:n_bootstrap) {
    # Bootstrap sample of species
    boot_indices <- sample(1:n_species, size = n_species, replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    # Impute bootstrap sample
    boot_imputed <- missForest(
      xmis = boot_data,
      verbose = FALSE
    )$ximp[, traitsSelect]
    
    bootstrap_results[[b]] <- boot_imputed
  }
  
  return(bootstrap_results)
}

# Apply bootstrap imputation
bootstrap_imputations <- bootstrap_imputation(traitsAux)
```

#### Option 3: Propagate Imputation Uncertainty Through Analysis

```r
# Function to run complete analysis pipeline on each imputed dataset
run_analysis_with_imputation <- function(imputed_datasets) {
  
  results_list <- list()
  
  for(i in 1:length(imputed_datasets)) {
    cat("Analyzing imputation", i, "\n")
    
    # Use i-th imputed dataset
    current_traits <- imputed_datasets[[i]]
    
    # Run PCA
    pca_result <- psych::principal(
      scale(current_traits), 
      nfactors = 4, 
      rotate = "varimax", 
      covar = TRUE
    )
    
    # Calculate TPDs
    trait_scores <- data.frame(pca_result$scores[, 1:4])
    sdTraits <- sqrt(diag(Hpi.diag(trait_scores)))
    
    tpds <- TPDsMean_large(
      species = rownames(trait_scores),
      means = trait_scores,
      sds = matrix(rep(sdTraits, nrow(trait_scores)), 
                   byrow = TRUE, ncol = 4),
      n_divisions = 30
    )
    
    # Store results
    results_list[[i]] <- list(
      pca = pca_result,
      tpds = tpds,
      variance_explained = pca_result$Vaccounted[2, ],
      loadings = pca_result$loadings
    )
  }
  
  return(results_list)
}

# Pool results and calculate uncertainty
pool_results <- function(analysis_results) {
  
  n_imputations <- length(analysis_results)
  
  # Pool variance explained
  variance_matrix <- sapply(analysis_results, function(x) x$variance_explained)
  pooled_variance <- rowMeans(variance_matrix)
  variance_uncertainty <- apply(variance_matrix, 1, sd)
  
  # Pool loadings
  loadings_array <- array(
    unlist(lapply(analysis_results, function(x) as.matrix(x$loadings))),
    dim = c(nrow(analysis_results[[1]]$loadings), 
            ncol(analysis_results[[1]]$loadings), 
            n_imputations)
  )
  
  pooled_loadings <- apply(loadings_array, c(1,2), mean)
  loadings_uncertainty <- apply(loadings_array, c(1,2), sd)
  
  return(list(
    variance_explained = pooled_variance,
    variance_se = variance_uncertainty,
    loadings = pooled_loadings,
    loadings_se = loadings_uncertainty,
    n_imputations = n_imputations
  ))
}
```

#### Option 4: Sensitivity Analysis Framework

```r
# Test sensitivity to imputation method
test_imputation_sensitivity <- function(trait_data) {
  
  methods <- list(
    missforest_default = function(x) missForest(x, verbose=FALSE)$ximp,
    missforest_conservative = function(x) missForest(x, ntree=500, maxiter=15, verbose=FALSE)$ximp,
    mice_pmm = function(x) complete(mice(x, method="pmm", m=1, printFlag=FALSE)),
    mice_cart = function(x) complete(mice(x, method="cart", m=1, printFlag=FALSE))
  )
  
  sensitivity_results <- list()
  
  for(method_name in names(methods)) {
    cat("Testing method:", method_name, "\n")
    
    imputed_data <- methods[[method_name]](trait_data)
    
    # Run key analyses
    pca_result <- psych::principal(scale(imputed_data[, traitsSelect]), 
                                   nfactors=4, rotate="varimax")
    
    sensitivity_results[[method_name]] <- list(
      variance_explained = pca_result$Vaccounted[2,],
      loadings = pca_result$loadings,
      method = method_name
    )
  }
  
  return(sensitivity_results)
}

# Compare methods
sensitivity_analysis <- test_imputation_sensitivity(traitsAux)
```

### Recommended Implementation Strategy

1. **Immediate fix**: Use multiple imputation (20+ datasets) with `missForest`
2. **Pool estimates**: Calculate means and standard errors across imputations  
3. **Propagate uncertainty**: Run key analyses on each imputed dataset
4. **Report uncertainty**: Include confidence intervals for all major findings
5. **Sensitivity testing**: Compare results across different imputation methods

### Impact on Study Conclusions

Proper uncertainty propagation may:
- Widen confidence intervals for variance explained by each component
- Affect statistical significance of PERMANOVA and dissimilarity analyses  
- Change interpretation of functional redundancy patterns
- Provide more honest assessment of analytical precision

---

## Computational Requirements

- Total runtime: ~45+ hours (dominated by null model script 05)
- Memory intensive: Large TPD calculations for multi-dimensional spaces
- Package dependencies: TPD, vegan, psych, missForest, V.PhyloMaker, and others