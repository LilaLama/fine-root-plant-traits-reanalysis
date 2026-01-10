# ----- Compare PCA loadings across imputation methods -----
# Extracts varimax-rotated loadings for all 10 traits across PC1-PC4
# Compares:
#   1. Original imputation (from script 06)
#   2. Mean imputation baseline (script 062) - naive approach (fill NAs with trait means)
#   3. Bootstrap imputation (1x OOB, script 061) - missForest with OOB-based uncertainty
#   4. Bootstrap imputation (2x OOB, script 063) - sensitivity analysis (doubled OOB noise)
#
# OUTPUT: data/Loadings_all_PCs_comparison.csv
# A table where each row is a trait (ph, ssd, sm, la, ln, sla, SRL, D, RTD, N)
# and each column is a loading value (correlation between original trait and principal component).
# Higher absolute loading = stronger contribution of that trait to that PC.
# Comparing loadings across methods shows whether imputation method affects PCA structure.
# Prerequisite: scripts 01, 061, 062, 063 have been run

# Load the four PCATotal objects
PCA_original <- readRDS("data/PCATotal_ImputedObs.rds") 
PCA_mean <- readRDS("data/PCATotal_mean_imputation.rds")  # Mean imputation baseline (script 062)
PCA_single_1x <- readRDS("data/imputed_bootstrap/PCATotal_ImputedObs_single.rds")  # Single run from 1x bootstrap (script 061)
PCA_single_2x <- readRDS("data/imputed_bootstrap_2x/PCATotal_ImputedObs_single.rds")  # Single run from 2x bootstrap (script 063)

# Check if objects loaded correctly
cat("PCA_original loaded:", !is.null(PCA_original), "\n")
cat("PCA_mean loaded:", !is.null(PCA_mean), "\n")
cat("PCA_single_1x loaded:", !is.null(PCA_single_1x), "\n")
cat("PCA_single_2x loaded:", !is.null(PCA_single_2x), "\n")

# Extract varimax-rotated loadings (10 traits × 4 PCs)
loadings_original <- as.matrix(PCA_original$PCA$loadings[, 1:4])
loadings_mean <- as.matrix(PCA_mean$PCA$loadings[, 1:4])
loadings_1x <- as.matrix(PCA_single_1x$PCA$loadings[, 1:4])
loadings_2x <- as.matrix(PCA_single_2x$PCA$loadings[, 1:4])

cat("\nLoadings_original dimensions:", dim(loadings_original), "\n")
cat("Loadings_original rownames:", paste(rownames(loadings_original), collapse=", "), "\n")
cat("Loadings_mean dimensions:", dim(loadings_mean), "\n")
cat("Loadings_1x dimensions:", dim(loadings_1x), "\n")
cat("Loadings_2x dimensions:", dim(loadings_2x), "\n")

cat("\n=== PCA Loadings Comparison: Original vs. Mean vs. 1x OOB vs. 2x OOB ===\n")

# Define trait order (as specified)
trait_order <- c("ph", "ssd", "sm", "la", "ln", "sla", "SRL", "D", "RTD", "N")

# Reorder loadings matrices by trait
loadings_original <- loadings_original[trait_order, ]
loadings_mean <- loadings_mean[trait_order, ]
loadings_1x <- loadings_1x[trait_order, ]
loadings_2x <- loadings_2x[trait_order, ]

# Create wide format comparison table with specified structure:
# Trait | PC1_Original | PC2_Original | PC3_Original | PC4_Original | 
#        | PC1_Mean | PC2_Mean | PC3_Mean | PC4_Mean |
#        | PC1_Boot1x | PC2_Boot1x | PC3_Boot1x | PC4_Boot1x | 
#        | PC1_Boot2x | PC2_Boot2x | PC3_Boot2x | PC4_Boot2x
comparison_wide <- data.frame(
  Trait = trait_order
)

# Add Original columns
for(pc in 1:4) {
  comparison_wide[[paste0("PC", pc, "_Original")]] <- round(loadings_original[, pc], 4)
}

# Add Mean imputation columns
for(pc in 1:4) {
  comparison_wide[[paste0("PC", pc, "_Mean")]] <- round(loadings_mean[, pc], 4)
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

# ===== VISUALIZATION SECTION =====
cat("\n=== Generating visualizations ===\n")

library(pheatmap)
library(RColorBrewer)

# === 1. BIPLOTS for each method ===
# Compute global limits for all biplots (same for all methods)
all_scores <- list(
  PCA_original$PCA$scores[, 1:4],
  PCA_mean$PCA$scores[, 1:4],
  PCA_single_1x$PCA$scores[, 1:4],
  PCA_single_2x$PCA$scores[, 1:4]
)

# Calculate global xlim and ylim for each component pair
comp_pairs <- list(c(1, 2), c(1, 3), c(3, 4), c(2, 3))
global_limits <- list()

for(k in seq_along(comp_pairs)) {
  c1 <- comp_pairs[[k]][1]
  c2 <- comp_pairs[[k]][2]
  
  # Get ranges across all four methods
  all_x <- sapply(all_scores, function(x) x[, c1])
  all_y <- sapply(all_scores, function(x) x[, c2])
  
  x_range <- range(all_x, na.rm = TRUE)
  y_range <- range(all_y, na.rm = TRUE)
  
  # Add 10% padding
  x_padding <- (x_range[2] - x_range[1]) * 0.1
  y_padding <- (y_range[2] - y_range[1]) * 0.1
  
  global_limits[[k]] <- list(
    xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
    ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)
  )
}

# Function to create biplot for each imputation method
create_biplot <- function(PCA_obj, main_title, outfile = NULL, global_limits) {
  loadings <- PCA_obj$PCA$loadings[, 1:4]
  scores <- PCA_obj$PCA$scores[, 1:4]
  
  # Create 2x2 grid: PC1 vs PC2, PC1 vs PC3, PC3 vs PC4, PC2 vs PC3
  png(outfile, width = 1400, height = 1200, res = 150)
  op <- par(mfrow = c(2, 2), mar = c(5, 5, 3, 2), oma = c(0, 0, 3, 0))
  
  comp_pairs <- list(c(1, 2), c(1, 3), c(3, 4), c(2, 3))
  pair_names <- c("PC1 vs PC2", "PC1 vs PC3", "PC3 vs PC4", "PC2 vs PC3")
  
  for(k in seq_along(comp_pairs)) {
    c1 <- comp_pairs[[k]][1]
    c2 <- comp_pairs[[k]][2]
    
    # Plot species scores with global limits
    plot(scores[, c1], scores[, c2],
         pch = 16, col = rgb(0, 0, 0, 0.3),
         xlab = paste0("PC", c1, " (", round(PCA_obj$Variance[c1] * 100, 1), "%)"),
         ylab = paste0("PC", c2, " (", round(PCA_obj$Variance[c2] * 100, 1), "%)"),
         main = pair_names[k],
         xlim = global_limits[[k]]$xlim,
         ylim = global_limits[[k]]$ylim,
         cex = 1.2)
    
    # Overlay trait loadings as arrows (scaled for visibility)
    arrow_scale <- 2.5
    arrows(0, 0, 
           loadings[, c1] * arrow_scale, loadings[, c2] * arrow_scale,
           col = rgb(1, 0, 0, 0.7), lwd = 2.5, length = 0.15)
    
    # Add trait labels
    text(loadings[, c1] * arrow_scale * 1.15, loadings[, c2] * arrow_scale * 1.15,
         rownames(loadings), col = rgb(1, 0, 0, 0.8), font = 2, cex = 1.1)
    
    # Add grid
    grid(col = rgb(0, 0, 0, 0.1))
  }
  
  # Add overall title
  mtext(main_title, outer = TRUE, cex = 1.5, font = 2)
  
  par(op)
  dev.off()
  cat("Saved:", outfile, "\n")
}

# Create biplots for each method
create_biplot(PCA_original, "Original missForest imputation", "Figures/Biplot_Original.png", global_limits)
create_biplot(PCA_mean, "Mean imputation", "Figures/Biplot_MeanImputation.png", global_limits)
create_biplot(PCA_single_1x, "Bootstrap 1x OOB", "Figures/Biplot_Bootstrap1xOOB.png", global_limits)
create_biplot(PCA_single_2x, "Bootstrap 2x OOB", "Figures/Biplot_Bootstrap2xOOB.png", global_limits)

# === 2. HEATMAP: Loadings comparison across methods ===
# Create a matrix with all loadings for each PC
for(pc in 1:4) {
  heatmap_data <- data.frame(
    Original = loadings_original[, pc],
    Mean = loadings_mean[, pc],
    Boot_1x = loadings_1x[, pc],
    Boot_2x = loadings_2x[, pc],
    row.names = trait_order
  )
  
  # Convert to matrix
  heatmap_mat <- as.matrix(heatmap_data)
  
  # Create heatmap
  png(paste0("Figures/Heatmap_Loadings_PC", pc, ".png"), width = 800, height = 600, res = 150)
  pheatmap(heatmap_mat,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           breaks = seq(-1, 1, length.out = 51),
           main = paste0("PC", pc, " Loadings Comparison across Imputation Methods"),
           display_numbers = TRUE,
           number_format = "%.3f",
           cellwidth = 60,
           cellheight = 25,
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  dev.off()
  cat("Saved: Figures/Heatmap_Loadings_PC", pc, ".png\n", sep = "")
}

# === 3. COMBINED HEATMAP: All PCs, all methods (trait × method×PC)
# Create a larger comparison matrix
combined_data <- data.frame()
for(pc in 1:4) {
  col_names <- paste0(c("Orig", "Mean", "1x", "2x"), "_PC", pc)
  combined_data <- cbind(
    combined_data,
    Original = loadings_original[, pc],
    Mean = loadings_mean[, pc],
    Boot_1x = loadings_1x[, pc],
    Boot_2x = loadings_2x[, pc]
  )
}
colnames(combined_data) <- apply(expand.grid(c("Original", "Mean", "1x OOB", "2x OOB"), 
                                              paste0("PC", 1:4)), 1, paste, collapse = "\n")
rownames(combined_data) <- trait_order

png("Figures/Heatmap_All_Loadings_Comparison.png", width = 1400, height = 700, res = 150)
pheatmap(as.matrix(combined_data),
         color = colorRampPalette(c("blue", "white", "red"))(50),
         breaks = seq(-1, 1, length.out = 51),
         main = "All PCA Loadings: Comparison across Imputation Methods & Components",
         display_numbers = TRUE,
         number_format = "%.2f",
         cellwidth = 35,
         cellheight = 25,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 9)
dev.off()
cat("Saved: Figures/Heatmap_All_Loadings_Comparison.png\n")

# === 4. SCATTER plots: Method comparison (e.g., Original vs. Mean)
png("Figures/Loadings_Scatter_OriginalVsMean.png", width = 800, height = 800, res = 150)
op <- par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
for(pc in 1:4) {
  plot(loadings_original[, pc], loadings_mean[, pc],
       pch = 16, col = rgb(0, 0, 1, 0.6), cex = 2,
       xlab = "Original (missForest)",
       ylab = "Mean imputation",
       main = paste0("PC", pc),
       xlim = c(-1, 1), ylim = c(-1, 1))
  text(loadings_original[, pc], loadings_mean[, pc],
       trait_order, pos = 3, cex = 0.9)
  # Add diagonal line (perfect agreement)
  abline(0, 1, col = "gray", lty = 2, lwd = 2)
  grid(col = rgb(0, 0, 0, 0.1))
}
par(op)
dev.off()
cat("Saved: Figures/Loadings_Scatter_OriginalVsMean.png\n")

# === 5. SCATTER plots: 1x OOB vs 2x OOB
png("Figures/Loadings_Scatter_1xVs2x.png", width = 800, height = 800, res = 150)
op <- par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
for(pc in 1:4) {
  plot(loadings_1x[, pc], loadings_2x[, pc],
       pch = 16, col = rgb(1, 0, 0, 0.6), cex = 2,
       xlab = "1x OOB",
       ylab = "2x OOB",
       main = paste0("PC", pc),
       xlim = c(-1, 1), ylim = c(-1, 1))
  text(loadings_1x[, pc], loadings_2x[, pc],
       trait_order, pos = 3, cex = 0.9)
  # Add diagonal line (perfect agreement)
  abline(0, 1, col = "gray", lty = 2, lwd = 2)
  grid(col = rgb(0, 0, 0, 0.1))
}
par(op)
dev.off()
cat("Saved: Figures/Loadings_Scatter_1xVs2x.png\n")

cat("\n=== All visualizations complete ===\n")
cat("Generated files in Figures/:\n")
cat("  - Biplot_Original.png, Biplot_MeanImputation.png, Biplot_Bootstrap1xOOB.png, Biplot_Bootstrap2xOOB.png\n")
cat("  - Heatmap_Loadings_PC1-4.png\n")
cat("  - Heatmap_All_Loadings_Comparison.png\n")
cat("  - Loadings_Scatter_OriginalVsMean.png\n")
cat("  - Loadings_Scatter_1xVs2x.png\n")

