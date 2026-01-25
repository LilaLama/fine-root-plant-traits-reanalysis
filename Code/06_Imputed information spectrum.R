### In this script we prepare plant trait data for analyses: 
# 0. Loading trait data (average of trait values at the species level from GROOT (fine root traits), which are already log-transformed, and loading of aboveground trait data from TRY (needs log-transformation). Taxonomical information for species with aboveground traits comes from The Plant List (package 'Taxonstand').
# 1. Selecting species with information for at least 50% of aboveground and 50% of fine-root traits. Acquisition of phylogeny for those species and imputation of missing trait values
# 2. Creating functional space based on those species (imputed + complete = "imputed dataset")
cat(paste0("\n\n Starting script #6 \n\n"))

########################################################################
### 0. Loading trait and taxonomy data at species level
########################################################################
rootTraits <- read.table("data/Root_traits.txt")[, c("SRL", "D", "RTD", "N")]
aboveTraits <- read.table("data/Above_traits.txt")
aboveTraits <- log10(aboveTraits)
aboveTaxonomy <- read.table("data/Above_taxonomy.txt")

########################################################################
### 1. Imputation of trait values
########################################################################
########################################################################
aboveInRoots <- intersect(rownames(aboveTraits), rownames(rootTraits))
plantTraitsAbove <- aboveTraits[aboveInRoots, ]
plantTraitsRoots <- rootTraits[aboveInRoots, ][rownames(plantTraitsAbove), ]
identical(rownames(plantTraitsRoots), rownames(plantTraitsAbove))
taxonomy <- aboveTaxonomy[aboveInRoots,]
AllTraitsAllInfo <- cbind(plantTraitsAbove, plantTraitsRoots, taxonomy) #1719 species
traitsSelect <- c("la", "ln", "ph", "sla", "ssd", "sm", "SRL", "D", "RTD", "N")
##### Lets select species with at least two traits below and 3 traits above measured (50% of traits):
traitsAbove <- c("la", "ln", "ph", "sla", "ssd", "sm")
traitsBelow <- c("SRL", "D", "RTD", "N")
maxNAAbove <- 3
maxNABelow <- 2
spSelectAbove <- which(rowSums(is.na(AllTraitsAllInfo[,traitsAbove])) <= maxNAAbove)
spSelectBelow <- which(rowSums(is.na(AllTraitsAllInfo[,traitsBelow])) <= maxNABelow)
spSelectALL <- spSelectBelow[which(names(spSelectBelow) %in% names(spSelectAbove))]
AllTraitsAllInfo <- AllTraitsAllInfo[names(spSelectALL), ] #1218 species
traitsUse <- AllTraitsAllInfo[, traitsSelect]

sp.list <- data.frame(species = gsub(pattern = "_", replacement = " ", 
                                     x = rownames(AllTraitsAllInfo)),
                      genus = AllTraitsAllInfo$genus,
                      family = AllTraitsAllInfo$family)
phyloPlants <- phylo.maker(sp.list = sp.list, tree = GBOTB.extended,
                           nodes = nodes.info.1,
                           output.sp.list = TRUE, output.tree = FALSE,
                           scenarios = "S1")
phylogenyAux <- phyloPlants$scenario.1
traitsUse <- traitsUse[phylogenyAux$tip.label, ]
identical(phylogenyAux$tip.label, rownames(traitsUse)) #checking all is ordered correctly
phylDissAux<- sqrt(cophenetic(phylogenyAux))
pcoaPhyl <- cmdscale(phylDissAux, k=10) 
traitsAux <- cbind(traitsUse, pcoaPhyl)
colnames(traitsAux) <- c(colnames(traitsUse), paste0("PC.", 1:ncol(pcoaPhyl)))
traitsAux <- traitsAux[rownames(AllTraitsAllInfo), ]
imputedTraits <- (missForest(xmis= traitsAux)$ximp[, traitsSelect])
identical(rownames(AllTraitsAllInfo), rownames(imputedTraits)) #checking all is ordered correctly
AllTraitsAllInfo[, traitsSelect] <- imputedTraits # Fill trait info with imputed values
saveRDS(AllTraitsAllInfo, "data/imputedTraits.rds")
saveRDS(traitsAux, "data/traitsAux.rds")      # For use in scripts 061, 063
saveRDS(traitsSelect, "data/traitsSelect.rds") # For use in scripts 061, 062, 063
saveRDS(traitsUse, "data/traitsUse.rds")      # For use in script 062

########################################################################
### 2. Functional space using imputed dataset
########################################################################
AllTraitsAllInfo <- readRDS("data/imputedTraits.rds")
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
#### Make sure orientation of all vectors is as in the paper:
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
                                        n_divisions = gridSize)
#Comp1. Comp3
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp3)))
PCATotal$TPDs2D$Comp1_Comp3 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp3), 
                                        means = PCATotal$traitsUse2D$Comp1_Comp3, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp3)),
                                                     byrow=T, ncol=2),
                                        n_divisions = gridSize)
#Comp1. Comp4
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp4)))
PCATotal$TPDs2D$Comp1_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp1_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp4)),
                                                     byrow=T, ncol=2),
                                        n_divisions = gridSize)

#Comp2. Comp3
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp2_Comp3)))
PCATotal$TPDs2D$Comp2_Comp3 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp2_Comp3), 
                                        means = PCATotal$traitsUse2D$Comp2_Comp3, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp2_Comp3)),
                                                     byrow=T, ncol=2),
                                        n_divisions = gridSize)
#Comp2. Comp4
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp2_Comp4)))
PCATotal$TPDs2D$Comp2_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp2_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp2_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp2_Comp4)),
                                                     byrow=T, ncol=2),
                                        n_divisions = gridSize)

#Comp3. Comp4
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp3_Comp4)))
PCATotal$TPDs2D$Comp3_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp3_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp3_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp3_Comp4)),
                                                     byrow=T, ncol=2),
                                        n_divisions = gridSize)
PCATotal$Readme <- "This object contains the functional space created by: 1. PCA (function psych::principal) using 1218 species with imputed information followed by varimax rotation. 2. PCA using the same set of species but without rotation. 3 Estimation of TPDs for 4 dimensions and between pairs of rotated components."


saveRDS(PCATotal, paste0("data/PCATotal_ImputedObs.rds"))

