### In this script we prepare plant trait data for analyses: 
# 0. Loading trait data (average of trait values at the species level from GROOT (fine root traits), which are already log-transformed, and loading of aboveground trait data from TRY (needs log-transformation). Taxonomical information for species with aboveground traits comes from The Plant List (package 'Taxonstand').
# 1. Selecting species with complete trait information for aboveground and fine-root traits separately, and for all traits combined.
# 2. Estimation of functional space (including estimation of dimensionality, estimation of PCA and varimax rotation, and estimation of TPD functions) based on complete information for aboveground traits, fine-root traits, and all traits combined.
cat(paste0("\n\n Starting script #1 \n\n"))
########################################################################
### 0. Loading trait and taxonomy data at species level
########################################################################
rootTraits <- read.table("data/Root_traits.txt")[, c("SRL", "D", "RTD", "N")]
aboveTraits <- read.table("data/Above_traits.txt")
aboveTraits <- log10(aboveTraits)
aboveTaxonomy <- read.table("data/Above_taxonomy.txt")

########################################################################
### 1. Selection of species with complete empirical information
########################################################################
rootTraitsComplete <- na.omit(rootTraits) # Species with complete information belowground (748)
aboveTraitsCompleteRows <- complete.cases(aboveTraits)
aboveTraitsComplete <- aboveTraits[aboveTraitsCompleteRows, ] # Species with complete information aboveground (2630)

taxonomy <- aboveTaxonomy[aboveTraitsCompleteRows, ]
aboveInRoots <- intersect(rownames(aboveTraitsComplete), rownames(rootTraitsComplete))
plantTraitsAbove <- aboveTraitsComplete[aboveInRoots, ]
plantTraitsRoots <- rootTraitsComplete[aboveInRoots, ][rownames(plantTraitsAbove), ]
identical(rownames(plantTraitsRoots), rownames(plantTraitsAbove))
AllTraitsAllInfo <- cbind(plantTraitsAbove, plantTraitsRoots) 
traitsUse <- AllTraitsAllInfo
AllTraitsAllInfoTax <- taxonomy[rownames(AllTraitsAllInfo),]
AllTraitsAllInfo <- cbind(traitsUse, AllTraitsAllInfoTax) # Species with complete information both above- and belowground (301)

########################################################################
### 2. Estimation of functional spaces for abovegroud, belowground, and combined
########################################################################

####################################
#### A. Only aboveground traits ####
gridSize <- 200
PCAAbove <- list()
PCAAbove$traits <- aboveTraitsComplete
PCAAbove$dimensions <- paran(PCAAbove$traits)$Retained
PCAAbove$means <- apply(PCAAbove$traits, 2, mean)
PCAAbove$sds <- apply(PCAAbove$traits, 2, sd)
PCAAbove$PCA <- psych::principal(scale(PCAAbove$traits), nfactors=PCAAbove$dimensions, 
                                 rotate="varimax", covar = T)
PCAAbove$Variance <- PCAAbove$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCAAbove$PCA$loadings**2))
for(i in 1:PCAAbove$dimensions){
  PCAAbove$PCA$scores[, i] <- PCAAbove$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCAAbove$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCAAbove$PCA$loadings))
PCAAbove$PCA$loadings <- PCAAbove$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCAAbove$traitsUse <- data.frame(PCAAbove$PCA$scores[, 1:PCAAbove$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCAAbove$traitsUse)))
PCAAbove$TPDs <- TPDsMean(species = rownames(PCAAbove$traitsUse), 
                    means = PCAAbove$traitsUse, 
                    sds = matrix(rep(sdTraits, nrow(PCAAbove$traitsUse)), byrow=T, 
                                 ncol=PCAAbove$dimensions),
                    n_divisions = gridSize)
saveRDS(PCAAbove, paste0("data/PCAAboveONLY_CompleteObs.rds"))
####################################
#### B. Only fine-root traits   ####
gridSize <- 200
PCABelow <- list()
plantTraitsBelow <- c("SRL", "D", "RTD", "N")
PCABelow$traits <- rootTraitsComplete
PCABelow$dimensions <- paran(PCABelow$traits)$Retained
PCABelow$means <- apply(PCABelow$traits, 2, mean)
PCABelow$sds <- apply(PCABelow$traits, 2, sd)
PCABelow$PCA <- psych::principal(scale(PCABelow$traits), nfactors=PCABelow$dimensions, 
                                 rotate="varimax", covar = T)
PCABelow$Variance <- PCABelow$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCABelow$PCA$loadings**2))
for(i in 1:PCABelow$dimensions){
  PCABelow$PCA$scores[, i] <- PCABelow$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCABelow$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCABelow$PCA$loadings))
PCABelow$PCA$loadings <- PCABelow$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCABelow$PCA$loadings**2)
PCABelow$traitsUse <- data.frame(PCABelow$PCA$scores[, 1:PCABelow$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCABelow$traitsUse)))
PCABelow$TPDs <- TPDsMean(species = rownames(PCABelow$traitsUse), 
                          means = PCABelow$traitsUse, 
                          sds = matrix(rep(sdTraits, nrow(PCABelow$traitsUse)), 
                                       byrow=T, ncol=PCABelow$dimensions),
                          n_divisions = gridSize)
saveRDS(PCABelow, paste0("data/PCABelowONLY_CompleteObs.rds")) 


##### PROCRUSTES ABOVE-BELOW:
commonSP <- intersect(rownames(PCABelow$traitsUse), rownames(PCAAbove$traitsUse))
vegan::protest(PCABelow$traitsUse[commonSP, ], PCAAbove$traitsUse[commonSP, ])

####################################
####     C. All traits         ##### 
AllTraits <- AllTraitsAllInfo[, c("la", "ln", "ph", "sla", "ssd", "sm",
                                  "SRL", "D", "RTD", "N")]
gridSize <- 30 # grid size (number of divisions per dimension, see TPD papers and package vignette) is smaller because of higher number of dimensions
PCATotal <- list()
PCATotal$traits <- AllTraits
PCATotal$dimensions <- paran(PCATotal$traits)$Retained
PCATotal$means <- apply(PCATotal$traits, 2, mean)
PCATotal$sds <- apply(PCATotal$traits, 2, sd)
PCATotal$AllInfo <- AllTraitsAllInfo
PCATotal$PCA <- psych::principal(scale(PCATotal$traits), nfactors=PCATotal$dimensions, 
                                 rotate="varimax", covar = T)
PCATotal$Variance <- PCATotal$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCATotal$PCA$loadings**2))
for(i in 1:PCATotal$dimensions){
  PCATotal$PCA$scores[, i] <- PCATotal$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCATotal$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCATotal$PCA$loadings))
PCATotal$PCA$loadings <- PCATotal$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCATotal$PCA$loadings**2)
#### Make sure orientation of all vectors is always consistent:
for(i in 1:PCATotal$dimensions){
  if(i == 1 & PCATotal$PCA$loadings["ph", i] < 0){ #ph is positive
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
  if(i == 2 & PCATotal$PCA$loadings["sla", i] < 0){ #sla is positive
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
  if(i == 3 & PCATotal$PCA$loadings["SRL", i] > 0){ #SRL is negative
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
  if(i == 4 & PCATotal$PCA$loadings["N", i] < 0){ #N is ppositive
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
}

#### UNROTATED PCA:
PCATotal$PCANoVarimax <- psych::principal(scale(PCATotal$traits), nfactors=PCATotal$dimensions,
                                          rotate="none", covar = T)
sqrtEigen2 <- sqrt(colSums(PCATotal$PCANoVarimax$loadings**2))
for(i in 1:PCATotal$dimensions){
  PCATotal$PCANoVarimax$scores[, i] <- PCATotal$PCANoVarimax$scores[, i] * sqrtEigen2[i] 
}
sqrtEigenMat <- matrix(rep(sqrtEigen2, nrow(PCATotal$PCANoVarimax$loadings)), byrow=T, 
                        nrow = nrow(PCATotal$PCANoVarimax$loadings))
PCATotal$PCANoVarimax$loadings <- PCATotal$PCANoVarimax$loadings / sqrtEigenMat
# NOTE: The results in PCATotal$PCANoVarimax are identical to those using princomp(scale(PCATotal$traits))
#### Make sure orientation of all vectors is as in the paper:
for(i in 1:PCATotal$dimensions){
  if(i == 1 & PCATotal$PCANoVarimax$loadings["ph", i] < 0){ #ph is positive
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
  }
  if(i == 2 & PCATotal$PCANoVarimax$loadings["sla", i] < 0){ #sla is positive
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
  }
  if(i == 3 & PCATotal$PCANoVarimax$loadings["SRL", i] > 0){ #SRL is negative
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
  }
  if(i == 4 & PCATotal$PCANoVarimax$loadings["N", i] < 0){ #N is ppositive
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 * PCATotal$PCANoVarimax$scores[, i]
  }
}

PCATotal$traitsUse <- data.frame(PCATotal$PCA$scores[, 1:PCATotal$dimensions]) 
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse)))
PCATotal$TPDs <- TPDsMean_large(species = rownames(PCATotal$traitsUse), 
                          means = PCATotal$traitsUse, 
                          sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse)), 
                                       byrow=T, ncol=PCATotal$dimensions),
                          n_divisions = gridSize) # TPD four dimensions

gridSize <- 200
PCATotal$traitsUse2D <- list()
PCATotal$TPDs2D <- list()
alphaUse <- 0.95
# TPD for the different planes combining pairs of dimensions
PCATotal$traitsUse2D$Comp1_Comp2 <- data.frame(PCATotal$PCA$scores[, 1:2]) 
PCATotal$traitsUse2D$Comp1_Comp3 <- data.frame(PCATotal$PCA$scores[, c(1, 3)]) 
PCATotal$traitsUse2D$Comp1_Comp4 <- data.frame(PCATotal$PCA$scores[, c(1, 4)]) 
PCATotal$traitsUse2D$Comp2_Comp3 <- data.frame(PCATotal$PCA$scores[, c(2, 3)]) 
PCATotal$traitsUse2D$Comp2_Comp4 <- data.frame(PCATotal$PCA$scores[, c(2, 4)]) 
PCATotal$traitsUse2D$Comp3_Comp4 <- data.frame(PCATotal$PCA$scores[, c(3, 4)]) 
#Comp1. Comp2
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp2)))
PCATotal$TPDs2D$Comp1_Comp2 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp2), 
                          means = PCATotal$traitsUse2D$Comp1_Comp2, 
                          sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp2)), 
                                       byrow=T, ncol=2),
                          alpha = alphaUse,
                          n_divisions = gridSize)
#Comp1. Comp3
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp3)))
PCATotal$TPDs2D$Comp1_Comp3 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp3), 
                                        means = PCATotal$traitsUse2D$Comp1_Comp3, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp3)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)
#Comp1. Comp4
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp4)))
PCATotal$TPDs2D$Comp1_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp1_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp4)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)

#Comp2. Comp3
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp2_Comp3)))
PCATotal$TPDs2D$Comp2_Comp3 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp2_Comp3), 
                                        means = PCATotal$traitsUse2D$Comp2_Comp3, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp2_Comp3)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)
#Comp2. Comp4
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp2_Comp4)))
PCATotal$TPDs2D$Comp2_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp2_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp2_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp2_Comp4)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)

#Comp3. Comp4
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp3_Comp4)))
PCATotal$TPDs2D$Comp3_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp3_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp3_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp3_Comp4)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)
PCATotal$Readme <- "This object contains the functional space created by: 1. PCA (function psych::principal) using 301 species with complete information followed by varimax rotation. 2. PCA using the same set of species but without rotation. 3 Estimation of TPDs for 4 dimensions and between pairs of rotated components."


saveRDS(PCATotal, paste0("data/PCATotal_CompleteObs.rds"))
