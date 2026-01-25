# ----- MEAN IMPUTATION: simple deterministic imputation for comparison -----
# Replace all missing trait values with the mean of observed values per trait.
# This is a naive baseline to compare against missForest imputation.
# Prerequisite: traitsUse, traitsSelect, AllTraitsAllInfo exist (from script 06)

library(psych)    # for principal() and paran()
library(ks)       # for Hpi.diag()

# Source auxiliary functions (needed for TPDs)
source("code/Aux_Functions.R")

# Load traitsUse, traitsSelect, and AllTraitsAllInfo from script 06
traitsUse <- readRDS("data/traitsUse.rds")
traitsSelect <- readRDS("data/traitsSelect.rds")
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")

# 1) Compute means per trait (only from observed values)
trait_means <- sapply(traitsUse[, traitsSelect], function(col){
  mean(col, na.rm = TRUE)
})

# 2) Impute: replace NAs with trait-wise means
imputedTraits_mean <- traitsUse[, traitsSelect]
for(j in colnames(imputedTraits_mean)){
  imputedTraits_mean[is.na(imputedTraits_mean[, j]), j] <- trait_means[j]
}

# 3) Create AllTraitsAllInfo with mean-imputed traits
AllTraitsAllInfo_mean <- AllTraitsAllInfo
AllTraitsAllInfo_mean[, traitsSelect] <- imputedTraits_mean
saveRDS(AllTraitsAllInfo_mean, file = "data/imputedTraits_mean.rds")

# 4) Run full PCA pipeline on mean-imputed data
# (copied from 06_Imputed information spectrum.R for consistency)

AllTraits_mean <- AllTraitsAllInfo_mean[, traitsSelect]
gridSize <- 30
PCATotal_mean <- list()
PCATotal_mean$traits <- AllTraits_mean
PCATotal_mean$dimensions <- paran(PCATotal_mean$traits)$Retained
PCATotal_mean$means <- apply(PCATotal_mean$traits, 2, mean)
PCATotal_mean$sds <- apply(PCATotal_mean$traits, 2, sd)
PCATotal_mean$AllInfo <- AllTraitsAllInfo_mean

PCATotal_mean$PCA <- psych::principal(scale(PCATotal_mean$traits), nfactors=PCATotal_mean$dimensions, 
                                      rotate="varimax", covar = T)
PCATotal_mean$Variance <- PCATotal_mean$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCATotal_mean$PCA$loadings**2))
for(i in 1:PCATotal_mean$dimensions){
  PCATotal_mean$PCA$scores[, i] <- PCATotal_mean$PCA$scores[, i] * sqrtEigen[i] 
}
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCATotal_mean$PCA$loadings)), byrow=TRUE, 
                       nrow = nrow(PCATotal_mean$PCA$loadings))
PCATotal_mean$PCA$loadings <- PCATotal_mean$PCA$loadings / sqrtEigenMat

# orientation (as in original script)
for(i in 1:PCATotal_mean$dimensions){
  if(i == 1 & PCATotal_mean$PCA$loadings["ph", i] < 0){
    PCATotal_mean$PCA$loadings[, i] <- -1 * PCATotal_mean$PCA$loadings[, i]
    PCATotal_mean$PCA$scores[, i] <- -1 *PCATotal_mean$PCA$scores[, i]
  }
  if(i == 2 & PCATotal_mean$PCA$loadings["sla", i] < 0){
    PCATotal_mean$PCA$loadings[, i] <- -1 * PCATotal_mean$PCA$loadings[, i]
    PCATotal_mean$PCA$scores[, i] <- -1 *PCATotal_mean$PCA$scores[, i]
  }
  if(i == 3 & PCATotal_mean$PCA$loadings["SRL", i] > 0){
    PCATotal_mean$PCA$loadings[, i] <- -1 * PCATotal_mean$PCA$loadings[, i]
    PCATotal_mean$PCA$scores[, i] <- -1 *PCATotal_mean$PCA$scores[, i]
  }
  if(i == 4 & PCATotal_mean$PCA$loadings["N", i] < 0){
    PCATotal_mean$PCA$loadings[, i] <- -1 * PCATotal_mean$PCA$loadings[, i]
    PCATotal_mean$PCA$scores[, i] <- -1 *PCATotal_mean$PCA$scores[, i]
  }
}

# unrotated
PCATotal_mean$PCANoVarimax <- psych::principal(scale(PCATotal_mean$traits), nfactors=PCATotal_mean$dimensions,
                                               rotate="none", covar = T)
sqrtEigen2 <- sqrt(colSums(PCATotal_mean$PCANoVarimax$loadings**2))
for(i in 1:PCATotal_mean$dimensions){
  PCATotal_mean$PCANoVarimax$scores[, i] <- PCATotal_mean$PCANoVarimax$scores[, i] * sqrtEigen2[i] 
}
sqrtEigenMat2 <- matrix(rep(sqrtEigen2, nrow(PCATotal_mean$PCANoVarimax$loadings)), byrow=TRUE, 
                        nrow = nrow(PCATotal_mean$PCANoVarimax$loadings))
PCATotal_mean$PCANoVarimax$loadings <- PCATotal_mean$PCANoVarimax$loadings / sqrtEigenMat2
for(i in 1:PCATotal_mean$dimensions){
  if(i == 1 & PCATotal_mean$PCANoVarimax$loadings["ph", i] < 0){
    PCATotal_mean$PCANoVarimax$loadings[, i] <- -1 * PCATotal_mean$PCANoVarimax$loadings[, i]
    PCATotal_mean$PCANoVarimax$scores[, i] <- -1 *PCATotal_mean$PCANoVarimax$scores[, i]
  }
  if(i == 2 & PCATotal_mean$PCANoVarimax$loadings["sla", i] < 0){
    PCATotal_mean$PCANoVarimax$loadings[, i] <- -1 * PCATotal_mean$PCANoVarimax$loadings[, i]
    PCATotal_mean$PCANoVarimax$scores[, i] <- -1 *PCATotal_mean$PCANoVarimax$scores[, i]
  }
  if(i == 3 & PCATotal_mean$PCANoVarimax$loadings["SRL", i] > 0){
    PCATotal_mean$PCANoVarimax$loadings[, i] <- -1 * PCATotal_mean$PCANoVarimax$loadings[, i]
    PCATotal_mean$PCANoVarimax$scores[, i] <- -1 *PCATotal_mean$PCANoVarimax$scores[, i]
  }
  if(i == 4 & PCATotal_mean$PCANoVarimax$loadings["N", i] < 0){
    PCATotal_mean$PCANoVarimax$loadings[, i] <- -1 * PCATotal_mean$PCANoVarimax$loadings[, i]
    PCATotal_mean$PCANoVarimax$scores[, i] <- -1 *PCATotal_mean$PCANoVarimax$scores[, i]
  }
}

# 4D TPDs
PCATotal_mean$traitsUse <- data.frame(PCATotal_mean$PCA$scores[, 1:PCATotal_mean$dimensions]) 
sdTraits <- sqrt(diag(Hpi.diag(PCATotal_mean$traitsUse)))
PCATotal_mean$TPDs <- TPDsMean_large(species = rownames(PCATotal_mean$traitsUse), 
                                     means = PCATotal_mean$traitsUse, 
                                     sds = matrix(rep(sdTraits, nrow(PCATotal_mean$traitsUse)), 
                                                  byrow=TRUE, ncol=PCATotal_mean$dimensions),
                                     n_divisions = gridSize)

# 2D TPDs
PCATotal_mean$traitsUse2D <- list()
PCATotal_mean$traitsUse2D$Comp1_Comp2 <- data.frame(PCATotal_mean$PCA$scores[, c(1,2)])
PCATotal_mean$traitsUse2D$Comp1_Comp3 <- data.frame(PCATotal_mean$PCA$scores[, c(1,3)])
PCATotal_mean$traitsUse2D$Comp1_Comp4 <- data.frame(PCATotal_mean$PCA$scores[, c(1,4)])
PCATotal_mean$traitsUse2D$Comp2_Comp3 <- data.frame(PCATotal_mean$PCA$scores[, c(2,3)])
PCATotal_mean$traitsUse2D$Comp2_Comp4 <- data.frame(PCATotal_mean$PCA$scores[, c(2,4)])
PCATotal_mean$traitsUse2D$Comp3_Comp4 <- data.frame(PCATotal_mean$PCA$scores[, c(3,4)])

PCATotal_mean$TPDs2D <- list()
PCATotal_mean$TPDs2D$Comp1_Comp2 <- TPDsMean(species = rownames(PCATotal_mean$traitsUse2D$Comp1_Comp2), 
                                             means = PCATotal_mean$traitsUse2D$Comp1_Comp2, 
                                             sds = matrix(c(sdTraits[1], sdTraits[2]), nrow(PCATotal_mean$traitsUse), 2, byrow = T),
                                             n_divisions = 50)
PCATotal_mean$TPDs2D$Comp1_Comp3 <- TPDsMean(species = rownames(PCATotal_mean$traitsUse2D$Comp1_Comp3), 
                                             means = PCATotal_mean$traitsUse2D$Comp1_Comp3, 
                                             sds = matrix(c(sdTraits[1], sdTraits[3]), nrow(PCATotal_mean$traitsUse), 2, byrow = T),
                                             n_divisions = 50)
PCATotal_mean$TPDs2D$Comp1_Comp4 <- TPDsMean(species = rownames(PCATotal_mean$traitsUse2D$Comp1_Comp4), 
                                             means = PCATotal_mean$traitsUse2D$Comp1_Comp4, 
                                             sds = matrix(c(sdTraits[1], sdTraits[4]), nrow(PCATotal_mean$traitsUse), 2, byrow = T),
                                             n_divisions = 50)
PCATotal_mean$TPDs2D$Comp2_Comp3 <- TPDsMean(species = rownames(PCATotal_mean$traitsUse2D$Comp2_Comp3), 
                                             means = PCATotal_mean$traitsUse2D$Comp2_Comp3, 
                                             sds = matrix(c(sdTraits[2], sdTraits[3]), nrow(PCATotal_mean$traitsUse), 2, byrow = T),
                                             n_divisions = 50)
PCATotal_mean$TPDs2D$Comp2_Comp4 <- TPDsMean(species = rownames(PCATotal_mean$traitsUse2D$Comp2_Comp4), 
                                             means = PCATotal_mean$traitsUse2D$Comp2_Comp4, 
                                             sds = matrix(c(sdTraits[2], sdTraits[4]), nrow(PCATotal_mean$traitsUse), 2, byrow = T),
                                             n_divisions = 50)
PCATotal_mean$TPDs2D$Comp3_Comp4 <- TPDsMean(species = rownames(PCATotal_mean$traitsUse2D$Comp3_Comp4), 
                                             means = PCATotal_mean$traitsUse2D$Comp3_Comp4, 
                                             sds = matrix(c(sdTraits[3], sdTraits[4]), nrow(PCATotal_mean$traitsUse), 2, byrow = T),
                                             n_divisions = 50)

PCATotal_mean$Readme <- "PCATotal from mean imputation: simple baseline for comparison."
saveRDS(PCATotal_mean, file = "data/PCATotal_mean_imputation.rds")

cat("Mean imputation done. Results saved to:\n")
cat("  - data/imputedTraits_mean.rds\n")
cat("  - data/PCATotal_mean_imputation.rds\n")
