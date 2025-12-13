# ----- Compare PCA loadings across imputation methods -----
# Extracts varimax-rotated loadings for all 10 traits across PC1-PC4
# Compares:
#   1. Single deterministic imputation (original complete cases)
#   2. Bootstrapped imputation (1x OOB)
#   3. Bootstrapped imputation (2x OOB)
# Prerequisite: scripts 01, 061, 063 have been run

# Load the three PCATotal objects (single deterministic runs from each imputation method)
PCA_original <- readRDS("data/PCATotal_ImputedObs.rds") 
PCA_single_1x <- readRDS("data/imputed_bootstrap/PCATotal_ImputedObs_single.rds")  # Single run from 1x bootstrap (script 061)
PCA_single_2x <- readRDS("data/imputed_bootstrap_2x/PCATotal_ImputedObs_single.rds")  # Single run from 2x bootstrap (script 063)

# Check if objects loaded correctly
cat("PCA_original loaded:", !is.null(PCA_original), "\n")
cat("PCA_single_1x loaded:", !is.null(PCA_single_1x), "\n")
cat("PCA_single_2x loaded:", !is.null(PCA_single_2x), "\n")

# Extract varimax-rotated loadings (10 traits × 4 PCs)
loadings_original <- as.matrix(PCA_original$PCA$loadings[, 1:4])
loadings_1x <- as.matrix(PCA_single_1x$PCA$loadings[, 1:4])
loadings_2x <- as.matrix(PCA_single_2x$PCA$loadings[, 1:4])

cat("\nLoadings_original dimensions:", dim(loadings_original), "\n")
cat("Loadings_original rownames:", paste(rownames(loadings_original), collapse=", "), "\n")
cat("Loadings_1x dimensions:", dim(loadings_1x), "\n")
cat("Loadings_2x dimensions:", dim(loadings_2x), "\n")

cat("\n=== PCA Loadings Comparison: Original vs. 1x OOB vs. 2x OOB ===\n")

# Define trait order (as specified)
trait_order <- c("ph", "ssd", "sm", "la", "ln", "sla", "SRL", "D", "RTD", "N")

# Reorder loadings matrices by trait
loadings_original <- loadings_original[trait_order, ]
loadings_1x <- loadings_1x[trait_order, ]
loadings_2x <- loadings_2x[trait_order, ]

# Create wide format comparison table with specified structure:
# Trait | PC1_Original | PC2_Original | PC3_Original | PC4_Original | PC1_Boot1x | PC2_Boot1x | PC3_Boot1x | PC4_Boot1x | PC1_Boot2x | PC2_Boot2x | PC3_Boot2x | PC4_Boot2x
comparison_wide <- data.frame(
  Trait = trait_order
)

# Add Original columns
for(pc in 1:4) {
  comparison_wide[[paste0("PC", pc, "_Original")]] <- round(loadings_original[, pc], 4)
}

# Add Bootstrap 1x OOB columns
for(pc in 1:4) {
  comparison_wide[[paste0("PC", pc, "_Boot1x")]] <- round(loadings_1x[, pc], 4)
}

# Add Bootstrap 2x OOB columns
for(pc in 1:4) {
  comparison_wide[[paste0("PC", pc, "_Boot2x")]] <- round(loadings_2x[, pc], 4)
}

# Print and save
cat("\n### Complete Loadings Comparison Table ###\n")
print(comparison_wide)

write.csv(comparison_wide, "data/Loadings_all_PCs_comparison.csv", row.names = FALSE)
saveRDS(comparison_wide, "data/Loadings_all_PCs_comparison.rds")

cat("\n=== Summary Statistics ===\n")
cat("Mean absolute difference (1x vs. original):", 
    round(mean(abs(loadings_1x - loadings_original)), 4), "\n")
cat("Mean absolute difference (2x vs. original):", 
    round(mean(abs(loadings_2x - loadings_original)), 4), "\n")
cat("Mean absolute difference (1x vs. 2x):", 
    round(mean(abs(loadings_1x - loadings_2x)), 4), "\n")

cat("\nComparison tables saved to data/Loadings_*.csv\n")
