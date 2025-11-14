### In this script we estimate what amount of functional space is occupied by each of teh considered groups (growth forms, families, biomes), and perform a null model to test whether this amount of space is smaller or larger than what would be expected by chance for the same number of species randomly selected, which is an indicator of the degree of functional redundancy among the species composing the group.
cat(paste0("\n\n Starting script #9 \n\n"))

start_time <- Sys.time()

data <- readRDS("data/PCATotal_ImputedObs.rds")
################################################################ 
################################################################ 
#### 1. Redundancy Woodiness #### 
woodiness <- read.table("data/woodiness_1218sp.txt", header = T)
dimensions <- data$dimensions
traitsUSE <- data$traitsUse
colnames(traitsUSE) <- paste0("Comp.",1:dimensions)
dataAnalysis <- traitsUSE
info <- data$AllInfo
nreps <- 10
classesWoody <- c("woody","non-woody")

TPDsAuxAbove <- data$TPDs2D$Comp1_Comp2
TPDsAuxBelow <- data$TPDs2D$Comp3_Comp4
TPDsAuxTotal <- data$TPDs
FRicWoodyAbove <- FRicWoodyBelow <- FRicWoodyTotal <- 
  matrix(NA, ncol = length(classesWoody), nrow = nreps + 1, 
         dimnames = list(c("Obs", paste0("Sim.", 1:nreps)), classesWoody))
commMatrixObs <- matrix(0, nrow = length(classesWoody), ncol = nrow(traitsUSE),
                     dimnames = list(classesWoody, rownames(traitsUSE)))
commMatrixObs["woody", which(woodiness == "woody")] <- 1
commMatrixObs["non-woody", which(woodiness == "non-woody")] <- 1
woodinessTPDc_Above_Obs <- TPDc(TPDs = TPDsAuxAbove, sampUnit = commMatrixObs)
FRicWoodyAbove[1, ] <-  TPDRichness(woodinessTPDc_Above_Obs)$communities$FRichness
woodinessTPDc_Below_Obs <- TPDc(TPDs = TPDsAuxBelow, sampUnit = commMatrixObs)
FRicWoodyBelow[1, ] <-  TPDRichness(woodinessTPDc_Below_Obs)$communities$FRichness
woodinessTPDc_Total_Obs <- TPDc_large(TPDs = TPDsAuxTotal, sampUnit = commMatrixObs)
FRicWoodyTotal[1, ] <-  TPDRichness_large(woodinessTPDc_Total_Obs)$communities$FRichness
for(i in 1:nreps){
  cat(paste0("\r REPETITION: ",i,"/",nreps,"\r"))
  ### Randomize species in each biome:
  commMatrixRand <- commMatrixObs
  commMatrixRand[,] <- 0
  for(GF in 1:nrow(commMatrixRand)){
    commMatrixRand[GF, sample(x = 1:ncol(commMatrixObs), size = rowSums(commMatrixObs)[GF])] <- 1
  }
  woodinessTPDc_Above_Rand <- TPDc(TPDs = TPDsAuxAbove, sampUnit = commMatrixRand)
  FRicWoodyAbove[1 + i, ] <-  TPDRichness(woodinessTPDc_Above_Rand)$communities$FRichness
  woodinessTPDc_Below_Rand <- TPDc(TPDs = TPDsAuxBelow, sampUnit = commMatrixRand)
  FRicWoodyBelow[1 + i, ] <-  TPDRichness(woodinessTPDc_Below_Rand)$communities$FRichness
  woodinessTPDc_Total_Rand <- TPDc_large(TPDs = TPDsAuxTotal, sampUnit = commMatrixRand)
  FRicWoodyTotal[1 + i, ] <-  TPDRichness_large(woodinessTPDc_Total_Rand)$communities$FRichness
}
saveRDS(FRicWoodyAbove, file = 'data/RedundancyWoodyAbove.rds')
saveRDS(FRicWoodyBelow, file = 'data/RedundancyWoodyBelow.rds')
saveRDS(FRicWoodyTotal, file = 'data/RedundancyWoodyTotal.rds')

################################################################ 
################################################################ 
#### 2. Redundancy families #### 
nspFam <- 15 #minimum number of species for the family to be considered
famPrint <- names(sort(table(data$AllInfo$family), decreasing = T)[which(sort(table(data$AllInfo$family), decreasing = T) >= nspFam)])
sum(table(data$AllInfo$family) >= nspFam)

FamMatTPD_Obs <- matrix(0, ncol = nrow(data$AllInfo), nrow = length(famPrint),
                        dimnames = list(famPrint, rownames(data$AllInfo)))
for(i in 1:nrow(FamMatTPD_Obs)){
  sp_i <- which(data$AllInfo$family == rownames(FamMatTPD_Obs)[i])
  FamMatTPD_Obs[i, sp_i] <- 1
}
####################################
# Null models Functional richness:
FRicFamiliesAbove <- FRicFamiliesBelow <- FRicFamiliesTotal <- matrix(NA, ncol = length(rownames(FamMatTPD_Obs)), nrow = nreps + 1,
                                              dimnames = list(c("Obs", paste0("Sim.", 1:nreps)),
                                                              rownames(FamMatTPD_Obs)))
FamsTPDc_Above_Obs <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = FamMatTPD_Obs)
FRicFamiliesAbove[1, ] <-  TPDRichness(FamsTPDc_Above_Obs)$communities$FRichness
FamsTPDc_Below_Obs <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = FamMatTPD_Obs)
FRicFamiliesBelow[1, ] <-  TPDRichness(FamsTPDc_Below_Obs)$communities$FRichness
FamsTPDc_Total_Obs <- TPDc_large(TPDs = data$TPDs, sampUnit = FamMatTPD_Obs)
FRicFamiliesTotal[1, ] <-  TPDRichness_large(FamsTPDc_Total_Obs)$communities$FRichness
for(i in 1:nreps){
  cat(paste0("\r FAMILIES; REPETITION: ",i,"/",nreps,"\r"))
  ### Randomize species in each Fam:
  FamMatTPD_Rand <- FamMatTPD_Obs
  FamMatTPD_Rand[,] <- 0
  for(Fam in 1:nrow(FamMatTPD_Rand)){
    FamMatTPD_Rand[Fam, sample(x = 1:ncol(FamMatTPD_Rand), size = rowSums(FamMatTPD_Obs)[Fam])] <- 1
  }
  FamsTPDc_Above_Rand <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = FamMatTPD_Rand)
  FRicFamiliesAbove[1 + i, ] <-  TPDRichness(FamsTPDc_Above_Rand)$communities$FRichness
  FamsTPDc_Below_Rand <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = FamMatTPD_Rand)
  FRicFamiliesBelow[1 + i, ] <-  TPDRichness(FamsTPDc_Below_Rand)$communities$FRichness
  FamsTPDc_Total_Rand <- TPDc_large(TPDs = data$TPDs, sampUnit = FamMatTPD_Rand)
  FRicFamiliesTotal[1 + i, ] <-  TPDRichness_large(FamsTPDc_Total_Rand)$communities$FRichness
}
saveRDS(FRicFamiliesAbove, file = 'data/RedundancyFamilyAbove.rds')
saveRDS(FRicFamiliesBelow, file = 'data/RedundancyFamilyBelow.rds')
saveRDS(FRicFamiliesTotal, file = 'data/RedundancyFamilyTotal.rds')

################################################################ 
################################################################ 
#### 2. Redundancy biomes #### 
#Biomes
biomesMat <- read.table("data/biomes_1218sp.txt", header = T)
nspBiom <- 15 #minimum number of species for the biome to be considered
BiomePrint <- names(which(colSums(biomesMat) >= nspBiom))
biomesMat <- biomesMat[, BiomePrint]
####################################
# Null models Functional richness:
FRicBiomesAbove <- FRicBiomesBelow <- FRicBiomesTotal <- matrix(NA, ncol = length(colnames(biomesMat)), nrow = nreps + 1,
                                              dimnames = list(c("Obs", paste0("Sim.", 1:nreps)),
                                                              colnames(biomesMat)))
bioMatTPD_Obs <- t(biomesMat)
biomesTPDc_Above_Obs <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = bioMatTPD_Obs)
FRicBiomesAbove[1, ] <-  TPDRichness(biomesTPDc_Above_Obs)$communities$FRichness
biomesTPDc_Below_Obs <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = bioMatTPD_Obs)
FRicBiomesBelow[1, ] <-  TPDRichness(biomesTPDc_Below_Obs)$communities$FRichness
biomesTPDc_Total_Obs <- TPDc_large(TPDs = data$TPDs, sampUnit = bioMatTPD_Obs)
FRicBiomesTotal[1, ] <-  TPDRichness_large(biomesTPDc_Total_Obs)$communities$FRichness
for(i in 1:nreps){
  cat(paste0("\r     BIOMES; REPETITION: ",i,"/",nreps,"\r"))
  ### Randomize species in each biome:
  bioMatTPD_Rand <- bioMatTPD_Obs
  bioMatTPD_Rand[,] <- 0
  for(biome in 1:nrow(bioMatTPD_Rand)){
    bioMatTPD_Rand[biome, sample(x = 1:ncol(bioMatTPD_Rand), size = rowSums(bioMatTPD_Obs)[biome])] <- 1
  }
  biomesTPDc_Above_Rand <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = bioMatTPD_Rand)
  FRicBiomesAbove[1 + i, ] <-  TPDRichness(biomesTPDc_Above_Rand)$communities$FRichness
  biomesTPDc_Below_Rand <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = bioMatTPD_Rand)
  FRicBiomesBelow[1 + i, ] <-  TPDRichness(biomesTPDc_Below_Rand)$communities$FRichness
  biomesTPDc_Total_Rand <- TPDc_large(TPDs = data$TPDs, sampUnit = bioMatTPD_Rand)
  FRicBiomesTotal[1 + i, ] <-  TPDRichness_large(biomesTPDc_Total_Rand)$communities$FRichness
}

saveRDS(FRicBiomesAbove, file = 'data/RedundancyBiomesAbove.rds')
saveRDS(FRicBiomesBelow, file = 'data/RedundancyBiomesBelow.rds')
saveRDS(FRicBiomesTotal, file = 'data/RedundancyBiomesTotal.rds')

end_time <- Sys.time()
cat(paste0("\n\n Script #9 completed in: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes \n\n"))