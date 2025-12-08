# ----- BOOTSTRAP / PARAMETRIC-BOOTSTRAP AROUND missForest IMPUTATION -----
# Prepare and run a single missForest imputation + report OOB errors per variable
# bootstrapped imputed datasets by adding parametric noise derived from OOB.
# The code tries to use objects produced by script 01 (aboveTraits, rootTraits).
# Adjust nboot / ntree / maxiter as needed.

# -------------------------------------------------------------------------
# OUTPUT FILES PRODUCED BY THIS SCRIPT:
# -------------------------------------------------------------------------
# data/imputedTraitsOOB.rds
#   - Complete missForest object (mf) with ximp, OOBerror, etc.
#
# data/imputed_bootstrap/imputed_single.rds
#   - AllTraitsAllInfo with the single (original) missForest imputation.
#
# data/imputed_bootstrap/imputed_boot_001.rds ... imputed_boot_050.rds
#   - 50 bootstrap versions of AllTraitsAllInfo, each with slightly perturbed
#     imputed values (parametric noise based on OOB errors).
#
# data/imputed_bootstrap/PCATotal_ImputedObs_single.rds
#   - PCATotal object for the single imputation (PCA scores, loadings, TPDs).
#
# data/imputed_bootstrap/PCATotal_ImputedObs_boot_001.rds ... _boot_050.rds
#   - PCATotal objects for each bootstrap replicate (PCA scores, loadings, TPDs).
#
# data/imputed_bootstrap/PCA_scores_boot_summary.rds
#   - List containing 'mean' and 'sd' of PCA scores across all bootstraps
#     (dimensions: species x 4 components).
# -------------------------------------------------------------------------


# ------------------ Replace simple imputation with parametric bootstrap propagation ------------------
# Prerequisite: traitsAux, traitsSelect, AllTraitsAllInfo exist as built above

library(missForest)

# 1) run single missForest imputation + report OOB errors per variable
mf <- missForest(xmis = traitsAux, variablewise = TRUE, verbose = TRUE)
ximp0 <- mf$ximp                      # full imputed matrix
imputedTraitsMatrix <- ximp0[, traitsSelect] # select real traits (excludes phylo)
names(mf$OOBerror) <- colnames(ximp0)
oob_raw <- mf$OOBerror                # per-variable OOB
saveRDS(mf, file = "data/imputedTraitsOOB.rds")

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

# 4) mark original missing positions (only add noise there)
missing_mat <- is.na(as.matrix(traitsAux)) # logical matrix of missing positions

# 5) Bootstrap settings
nboot <- 50
outdir <- "data/imputed_bootstrap"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# helper-function to create PCATotal from AllTraitsAllInfo object (copied from 06_Imputed information spectrum.R line 60ff)
'''
This helper function encapsulates the functional-space construction workflow from script 06_Imputed information spectrum.R 
so it can be reused on different imputed datasets (e.g., bootstrap replicates). 
It takes an AllTraitsAllInfo_obj (a data frame with imputed trait values plus taxonomic metadata) 
and an optional suffix for saving, and returns a PCATotal list containing PCA results, TPD estimates, and summary statistics.
'''
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
  PCATotal$Readme <- "PCATotal generated from imputed AllTraitsAllInfo"
  # save if suffix supplied
  if(!is.null(suffix)){
    saveRDS(PCATotal, file.path(outdir, paste0("PCATotal_ImputedObs", suffix, ".rds")))
  }
  return(PCATotal)
}

############ CONTINUE HERE ############

# 6) Create and save nboot imputations and PCATotal for each; also store PCA scores for uncertainty summaries
scores_list <- vector("list", nboot + 1) # slot 1 = single deterministic imputation, slots 2..nboot+1 = bootstraps
# boot 0: single imputation (ximp0)
imputed0 <- ximp0[, traitsSelect]
AllTraitsAllInfo_single <- AllTraitsAllInfo
AllTraitsAllInfo_single[, traitsSelect] <- imputed0
saveRDS(AllTraitsAllInfo_single, file = file.path(outdir, "imputed_single.rds"))

# Single imputation: compute PCATotal and extract scores
PCA_single <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_single, suffix = "_single")
scores_list[[1]] <- PCA_single$PCA$scores

# iterative boots
cat("\nBootstrapping imputations...\n")
pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
for(b in seq_len(nboot)){
  utils::setTxtProgressBar(pb, b); flush.console()
  set.seed(10000 + b)
  xboot <- ximp0
  for(j in seq_len(ncol(xboot))){
    miss_idx <- which(missing_mat[, j])
    if(length(miss_idx) > 0){
      xboot[miss_idx, j] <- rnorm(length(miss_idx), mean = ximp0[miss_idx, j], sd = noise_sd[j])
    }
  }
  imputed_b <- xboot[, traitsSelect]
  AllTraitsAllInfo_b <- AllTraitsAllInfo
  AllTraitsAllInfo_b[, traitsSelect] <- imputed_b
  saveRDS(AllTraitsAllInfo_b, file = file.path(outdir, paste0("imputed_boot_", sprintf("%03d", b), ".rds")))
  PCAb <- make_PCAtotal_from_AllTraits(AllTraitsAllInfo_b, suffix = paste0("_boot_", sprintf("%03d", b)))
  scores_list[[b + 1]] <- PCAb$PCA$scores
}
close(pb); cat("\nBootstrapping done.\n")

# 7) Aggregate PCA scores over bootstraps: mean and SD per species x component
scores_array <- simplify2array(lapply(scores_list, function(x) x[, 1:4])) # species x comps x boots
if(length(dim(scores_array)) == 3){
  score_mean <- apply(scores_array, c(1,2), mean)
  score_sd   <- apply(scores_array, c(1,2), sd)
  saveRDS(list(mean = score_mean, sd = score_sd), file = file.path(outdir, "PCA_scores_boot_summary.rds"))
}

# --- Compare single (deterministic) PCA scores vs. bootstrap PCs ---
if(length(scores_list) > 1){
  single_scores <- scores_list[[1]][, 1:4, drop = FALSE]
  boot_array    <- simplify2array(lapply(scores_list[-1], function(x) x[, 1:4]))
  if(length(dim(boot_array)) == 3){
    boot_mean <- apply(boot_array, c(1,2), mean)
    boot_sd   <- apply(boot_array, c(1,2), sd)
    delta     <- boot_mean - single_scores
    rel_ratio <- delta / (abs(single_scores) + 1e-8)  # relative change (avoid div/0)
    saveRDS(list(single = single_scores,
                 boot_mean = boot_mean,
                 boot_sd = boot_sd,
                 delta = delta,
                 rel_ratio = rel_ratio),
            file = file.path(outdir, "PCA_scores_single_vs_boot.rds"))
  }
}

# Visualization

