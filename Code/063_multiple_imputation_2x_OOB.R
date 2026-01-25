# ----- BOOTSTRAP WITH 2x OOB: sensitivity analysis -----
# Same as 063_multiple_imputation_incrOOB.R, but with noise_sd multiplied by 2
# to explore sensitivity to OOB-based uncertainty scaling.
# Results saved to data/imputed_bootstrap_2x/ to avoid overwriting original.
# Prerequisite: traitsAux, traitsSelect, AllTraitsAllInfo exist (from script 06)

library(missForest)
library(psych)    # for principal() and paran()
library(ks)       # for Hpi.diag()

# Source auxiliary functions (TPDsMean_large, etc.)
source("code/Aux_Functions.R")

# Load traitsAux and traitsSelect from script 06
traitsAux <- readRDS("data/traitsAux.rds")
traitsSelect <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")

# 1) run single missForest imputation + report OOB errors per variable
mf <- missForest(xmis = traitsAux, variablewise = TRUE, verbose = TRUE)
ximp0 <- mf$ximp                      # full imputed matrix
imputedTraitsMatrix <- ximp0[, traitsSelect] # select real traits (excludes phylo)
names(mf$OOBerror) <- colnames(ximp0)
oob_raw <- mf$OOBerror * 2 # per-variable OOB, scaled by 2.0 for sensitivity analysis
saveRDS(mf, file = "data/imputedTraits_2xOOB.rds")

# 2) extract OOB only for real traits (exclude phylo PCs)
oob_norm <- oob_raw[traitsSelect]

# 3) compute noise SD per variable
# sd_obs = observed SD per trait (log10-scale), using only non-NA values from the original data;
sd_obs <- sapply(colnames(ximp0), function(j){
  v <- traitsAux[, j][!is.na(traitsAux[, j])]
  if(length(v) > 1) sd(v) else 0
})
noise_sd <- oob_norm * sd_obs # still log10 scale
noise_sd[is.na(noise_sd) | noise_sd <= 0] <- 1e-8 # replace non-positive or missing values with a very small positive number

cat("\n*** OOB-BASED NOISE MULTIPLIED BY 2.0 ***\n")
cat("Noise SDs per trait (log10 scale):\n")
print(noise_sd[traitsSelect])

# 4) mark original missing positions (only add noise there)
missing_mat <- is.na(as.matrix(traitsAux)) # logical matrix of missing positions

# 5) Bootstrap settings
nboot <- 50
outdir <- "data/imputed_bootstrap_2x"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# helper-function to create PCATotal from AllTraitsAllInfo object (copied from 06_Imputed information spectrum.R line 60ff)
make_PCAtotal_from_AllTraits <- function(AllTraitsAllInfo_obj, suffix = NULL){
  AllTraits <- AllTraitsAllInfo_obj[, traitsSelect]
  gridSize <- 30
  PCATotal <- list()
  PCATotal$traits <- AllTraits
  PCATotal$dimensions <- paran(PCATotal$traits)$Retained
  PCATotal$means <- apply(PCATotal$traits, 2, mean)
  PCATotal$sds <- apply(PCATotal$traits, 2, sd)
  PCATotal$AllInfo <- AllTraitsAllInfo_obj
  PCATotal$PCA <- psych::principal(scale(PCATotal$traits), nfactors=PCATotal$dimensions, 
                                   rotate="varimax", covar = T)
  PCATotal$Variance <- PCATotal$PCA$Vaccounted[2,]
  sqrtEigen <- sqrt(colSums(PCATotal$PCA$loadings**2))
  for(i in 1:PCATotal$dimensions){
    PCATotal$PCA$scores[, i] <- PCATotal$PCA$scores[, i] * sqrtEigen[i] 
  }
  sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCATotal$PCA$loadings)), byrow=TRUE, 
                         nrow = nrow(PCATotal$PCA$loadings))
  PCATotal$PCA$loadings <- PCATotal$PCA$loadings / sqrtEigenMat
  # orientation (as in original script)
  for(i in 1:PCATotal$dimensions){
    if(i == 1 & PCATotal$PCA$loadings["ph", i] < 0){
      PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
      PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
    }
    if(i == 2 & PCATotal$PCA$loadings["sla", i] < 0){
      PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
      PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
    }
    if(i == 3 & PCATotal$PCA$loadings["SRL", i] > 0){
      PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
      PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
    }
    if(i == 4 & PCATotal$PCA$loadings["N", i] < 0){
      PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
      PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
    }
  }
  # unrotated
  PCATotal$PCANoVarimax <- psych::principal(scale(PCATotal$traits), nfactors=PCATotal$dimensions,
                                            rotate="none", covar = T)
  sqrtEigen2 <- sqrt(colSums(PCATotal$PCANoVarimax$loadings**2))
  for(i in 1:PCATotal$dimensions){
    PCATotal$PCANoVarimax$scores[, i] <- PCATotal$PCANoVarimax$scores[, i] * sqrtEigen2[i] 
  }
  sqrtEigenMat2 <- matrix(rep(sqrtEigen2, nrow(PCATotal$PCANoVarimax$loadings)), byrow=TRUE, 
                          nrow = nrow(PCATotal$PCANoVarimax$loadings))
  PCATotal$PCANoVarimax$loadings <- PCATotal$PCANoVarimax$loadings / sqrtEigenMat2
  for(i in 1:PCATotal$dimensions){
    if(i == 1 & PCATotal$PCANoVarimax$loadings["ph", i] < 0){
      PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
      PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
    }
    if(i == 2 & PCATotal$PCANoVarimax$loadings["sla", i] < 0){
      PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
      PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
    }
    if(i == 3 & PCATotal$PCANoVarimax$loadings["SRL", i] > 0){
      PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
      PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
    }
    if(i == 4 & PCATotal$PCANoVarimax$loadings["N", i] < 0){
      PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
      PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
    }
  }
  # 4D TPDs (as before)
  PCATotal$traitsUse <- data.frame(PCATotal$PCA$scores[, 1:PCATotal$dimensions]) 
  sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse)))
  PCATotal$TPDs <- TPDsMean_large(species = rownames(PCATotal$traitsUse), 
                                  means = PCATotal$traitsUse, 
                                  sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse)), 
                                               byrow=TRUE, ncol=PCATotal$dimensions),
                                  n_divisions = gridSize)
  PCATotal$Readme <- "PCATotal generated from imputed AllTraitsAllInfo (2x OOB noise)"
  # save if suffix supplied
  if(!is.null(suffix)){
    saveRDS(PCATotal, file.path(outdir, paste0("PCATotal_ImputedObs", suffix, ".rds")))
  }
  return(PCATotal)
}

# 6) Create and save nboot imputations and PCATotal for each; also store PCA scores and loadings for uncertainty summaries
scores_list <- vector("list", nboot + 1) # slot 1 = single deterministic imputation, slots 2..nboot+1 = bootstraps
loadings_list <- vector("list", nboot + 1) # store loadings for each bootstrap
# boot 0: single imputation (ximp0)
imputed0 <- ximp0[, traitsSelect]
AllTraitsAllInfo_single <- AllTraitsAllInfo
AllTraitsAllInfo_single[, traitsSelect] <- imputed0
saveRDS(AllTraitsAllInfo_single, file = file.path(outdir, "imputed_single.rds"))

# Single imputation: compute PCATotal and extract scores + loadings
PCA_single <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_single, suffix = "_single")
scores_list[[1]] <- PCA_single$PCA$scores
loadings_list[[1]] <- as.matrix(PCA_single$PCA$loadings[, 1:4])

# iterative boots
cat("\nBootstrapping imputations (2x OOB noise)...\n")
pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
for(b in seq_len(nboot)){
  utils::setTxtProgressBar(pb, b); flush.console()
  set.seed(10000 + b)
  xboot <- ximp0
  for(j in seq_len(ncol(xboot))){
    miss_idx <- which(missing_mat[, j])
    if(length(miss_idx) > 0){
      xboot[miss_idx, j] <- rnorm(length(miss_idx), mean = ximp0[miss_idx, j], sd = noise_sd[j]) # 2x amplified noise
    }
  }
  imputed_b <- xboot[, traitsSelect]
  AllTraitsAllInfo_b <- AllTraitsAllInfo
  AllTraitsAllInfo_b[, traitsSelect] <- imputed_b
  saveRDS(AllTraitsAllInfo_b, file = file.path(outdir, paste0("imputed_boot_", sprintf("%03d", b), ".rds")))
  PCAb <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_b, suffix = paste0("_boot_", sprintf("%03d", b)))
  scores_list[[b + 1]] <- PCAb$PCA$scores
  loadings_list[[b + 1]] <- as.matrix(PCAb$PCA$loadings[, 1:4])
}
close(pb); cat("\nBootstrapping done (2x OOB).\n")

# 7) Aggregate PCA scores over bootstraps: mean and SD per species x component
scores_array <- simplify2array(lapply(scores_list, function(x) x[, 1:4])) # species x comps x boots
if(length(dim(scores_array)) == 3){
  score_mean <- apply(scores_array, c(1,2), mean)
  score_sd   <- apply(scores_array, c(1,2), sd)
  saveRDS(list(mean = score_mean, sd = score_sd), file = file.path(outdir, "PCA_scores_boot_summary.rds"))
}

# 8) Aggregate PCA loadings over bootstraps: mean and SD per trait x component
loadings_array <- simplify2array(loadings_list) # traits x comps x boots
if(length(dim(loadings_array)) == 3){
  loadings_mean <- apply(loadings_array, c(1,2), mean)
  loadings_sd   <- apply(loadings_array, c(1,2), sd)
  colnames(loadings_mean) <- paste0("PC", 1:4)
  colnames(loadings_sd) <- paste0("PC", 1:4)
  saveRDS(list(mean = loadings_mean, sd = loadings_sd, all_boots = loadings_array), 
          file = file.path(outdir, "PCA_loadings_boot_summary.rds"))
  # Also save as CSV for easy inspection
  loadings_summary_df <- data.frame(
    Trait = rownames(loadings_mean),
    PC1_mean = loadings_mean[,1], PC1_sd = loadings_sd[,1],
    PC2_mean = loadings_mean[,2], PC2_sd = loadings_sd[,2],
    PC3_mean = loadings_mean[,3], PC3_sd = loadings_sd[,3],
    PC4_mean = loadings_mean[,4], PC4_sd = loadings_sd[,4]
  )
  write.csv(loadings_summary_df, file = file.path(outdir, "PCA_loadings_boot_summary.csv"), row.names = FALSE)
  cat("\nPCA loadings aggregated across", nboot + 1, "imputations (single + bootstraps).\n")
  cat("Saved to:", file.path(outdir, "PCA_loadings_boot_summary.rds"), "\n")
}

# --- Compare bootstrap PCs vs. bootstrap PCs (2x)  ---
# Load 1x OOB results for comparison
results_1x <- readRDS("data/imputed_bootstrap/PCA_scores_single_vs_boot.rds")
boot_mean_1x <- results_1x$boot_mean
boot_sd_1x <- results_1x$boot_sd

# 2x OOB results
if(length(scores_list) > 1){
  single_scores <- scores_list[[1]][, 1:4, drop = FALSE]
  boot_array    <- simplify2array(lapply(scores_list[-1], function(x) x[, 1:4]))
  if(length(dim(boot_array)) == 3){
    boot_mean_2x <- apply(boot_array, c(1,2), mean)
    boot_sd_2x   <- apply(boot_array, c(1,2), sd)
    delta     <- boot_mean_2x - boot_mean_1x
    rel_ratio <- delta / (abs(boot_mean_1x) + 1e-8)  # relative change (avoid div/0)
    saveRDS(list(boot_mean_1x = boot_mean_1x,
                 boot_sd_1x = boot_sd_1x,
                 boot_mean_2x = boot_mean_2x,
                 boot_sd_2x = boot_sd_2x,
                 delta = delta,
                 rel_ratio = rel_ratio),
            file = file.path(outdir, "PCA_scores_1x_vs_2x_boot.rds"))

    # Visualization: 1x OOB vs. 2x OOB bootstrap means with error bars
    comp_pairs <- list(c(1,2), c(3,4))
    pair_names <- c("PC1_vs_PC2", "PC3_vs_PC4")
    for(k in seq_along(comp_pairs)){
      c1 <- comp_pairs[[k]][1]; c2 <- comp_pairs[[k]][2]
      png(file.path(outdir, paste0("PCA_scores_1x_vs_2x_", pair_names[k], ".png")),
          width = 1200, height = 800, res = 150)
      op <- par(mfrow = c(1,1), mar = c(5,5,2,2))
      plot(boot_mean_1x[, c1], boot_mean_1x[, c2],
           pch = 16, col = rgb(0,0,0,0.5),
           xlab = paste0("Bootstrap mean (1x OOB) PC", c1),
           ylab = paste0("Bootstrap mean (1x OOB) PC", c2),
           main = paste0("1x OOB vs. 2x OOB bootstrap (", pair_names[k], ")"))
      points(boot_mean_2x[, c1], boot_mean_2x[, c2],
             pch = 16, col = rgb(1,0,0,0.5))
      # error bars ±1 SD for 2x
      segments(boot_mean_2x[, c1] - boot_sd_2x[, c1], boot_mean_2x[, c2],
               boot_mean_2x[, c1] + boot_sd_2x[, c1], boot_mean_2x[, c2],
               col = rgb(1,0,0,0.3))
      segments(boot_mean_2x[, c1], boot_mean_2x[, c2] - boot_sd_2x[, c2],
               boot_mean_2x[, c1], boot_mean_2x[, c2] + boot_sd_2x[, c2],
               col = rgb(1,0,0,0.3))
      legend("topright",
             legend = c("Bootstrap mean (1x OOB)", "Bootstrap mean ±1 SD (2x OOB)"),
             pch = c(16,16), col = c(rgb(0,0,0,0.5), rgb(1,0,0,0.5)), bty = "n")
      par(op)
      dev.off()
    }
  }
}

cat("\n2x OOB sensitivity analysis complete.\n")
cat("Results saved to: data/imputed_bootstrap_2x/\n")
cat("Compare with: data/imputed_bootstrap/ (original 1x OOB)\n")
