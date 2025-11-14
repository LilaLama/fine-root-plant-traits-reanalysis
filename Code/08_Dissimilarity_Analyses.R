### In this script we perform a dissimilarity analysis between pairs of growth forms (woody vs herbaceous species), biomes and families based on the occupation of the aboveground plane and belowground plane: 
cat(paste0("\n\n Starting script #8 \n\n"))


data <- readRDS("data/PCATotal_ImputedObs.rds")

################################################################ 
################################################################ 
#### 1. Dissimilarity between woody and herbaceous species #### 
woodiness <- read.table("data/woodiness_1218sp.txt", header = T)
nonwoodySP <- which(woodiness == "non-woody")
woodySP <- which(woodiness == "woody")
commNames <- c("woody", "nonwoody")
commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(data$traitsUse),
                     dimnames = list(commNames, rownames(data$traitsUse)))
commMatrix["woody", woodySP] <- 1
commMatrix["nonwoody", nonwoodySP] <- 1

TPDcAuxAbove <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = commMatrix)
dissimilaritiesWoodinessAbove <- dissim(TPDcAuxAbove)$communities$dissimilarity[1,2]
TPDcAuxBelow <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = commMatrix)
dissimilaritiesWoodinessBelow <- dissim(TPDcAuxBelow)$communities$dissimilarity[1,2]

################################################################ 
################################################################ 
#### 2. Dissimilarity between families #### 
nspFam <- 15 #minimum number of species for the family to be considered
famPrint <- names(sort(table(data$AllInfo$family), decreasing = T)[which(sort(table(data$AllInfo$family), decreasing = T) >= nspFam)])
sum(table(data$AllInfo$family) >= nspFam) # Families with at least nspFam species
famMat <- matrix(0, ncol = nrow(data$traitsUse), nrow = length(famPrint),
                 dimnames = list(famPrint, rownames(data$traitsUse)))

for(i in 1:length(famPrint)){
    sp_i <- which(data$AllInfo$family == famPrint[i])
    famMat[i, sp_i] <- 1
  }
TPDcFamAbove <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = famMat)
dissimFamilyAbove <- dissim(TPDcFamAbove)$communities$dissimilarity
TPDcFamBelow <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = famMat)
dissimFamilyBelow <- dissim(TPDcFamBelow)$communities$dissimilarity

saveRDS(dissimFamilyAbove, "data/dissim_Family_Above_imputedDataset.rds")
saveRDS(dissimFamilyBelow, "data/dissim_Family_Below_imputedDataset.rds")


################################################################ 
################################################################ 
#### 3. Dissimilarity between biomes #### 
#Biomes data
biomesMat <- read.table("data/biomes_1218sp.txt", header = T)
nspBiom <- 15 #minimum number of species for the biome to be considered
BiomePrint <- names(which(colSums(biomesMat) >= nspBiom))
BiomeMat <- t(biomesMat[, BiomePrint])

TPDcBiomeAbove <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = BiomeMat)
dissimBiomeAbove <- dissim(TPDcBiomeAbove)$communities$dissimilarity
TPDcBiomeBelow <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = BiomeMat)
dissimBiomeBelow <- dissim(TPDcBiomeBelow)$communities$dissimilarity

saveRDS(dissimBiomeAbove, "data/dissim_Biome_Above_imputedDataset.rds")
saveRDS(dissimBiomeBelow, "data/dissim_Biome_Below_imputedDataset.rds")


###### CLIMATE DISSIMILARITY BETWEEN BIOMES:
climBiomes <- as.data.frame(readRDS("data/Climate_biomes.rds"))
dissim_biome_Clmate <- as.matrix(FD::gowdis(climBiomes))
colnames(dissim_biome_Clmate)<-rownames(dissim_biome_Clmate)<-gsub(" ",".",rownames(dissim_biome_Clmate))
colnames(dissim_biome_Clmate)<-rownames(dissim_biome_Clmate)<-gsub("/",".",rownames(dissim_biome_Clmate))

dissim_biome_Clmate <- dissim_biome_Clmate[BiomePrint, BiomePrint]


###########################################################
#### Analyses ####
#### 
limX <- c(0, 1)
limY <- c(0, 1)
#mantel test families
par(mfrow=c(2,2))
vegan::mantel(dissimFamilyAbove, dissimFamilyBelow)
plot(as.dist(dissimFamilyBelow), as.dist(dissimFamilyAbove), pch = 21, 
     bg = rgb(.2, .3, .3, alpha = 0.5), col = NA, cex = 1.5,
     axes=T, main = "", xlim = limX, ylim = limY,
     xlab = "Dissimilarity fine-root plane",
     ylab = "Dissimilarity aboveground plane")
abline(0, 1, lty = 2, col = "grey60")
mod<-lmodel2(as.dist(dissimFamilyAbove)~as.dist(dissimFamilyBelow))
curve(expr= mod$regression.results$Intercept[2]+ mod$regression.results$Slope[2]*x,
      add=T, col=rgb(.2, .3, .3), lwd=5, lty=1)

#mantel test Biomes
vegan::mantel(dissimBiomeAbove, dissimBiomeBelow)
plot(as.dist(dissimBiomeBelow), as.dist(dissimBiomeAbove), pch = 21, 
     bg = rgb(.2, .3, .3, alpha = 0.5), col = NA, cex = 1.5,
     axes=T, main = "", xlim = limX, ylim = limY,
     xlab = "Dissimilarity fine-root plane",
     ylab = "Dissimilarity aboveground plane")
abline(0, 1, lty = 2, col = "grey60")
mod<-lmodel2(as.dist(dissimBiomeAbove)~as.dist(dissimBiomeBelow))
curve(expr= mod$regression.results$Intercept[2]+ mod$regression.results$Slope[2]*x,
      add=T, col=rgb(.2, .3, .3), lwd=5, lty=1)

#Climate
vegan::mantel(dissimBiomeAbove, dissim_biome_Clmate)
plot(as.dist(dissim_biome_Clmate), as.dist(dissimBiomeAbove), pch = 21, 
     bg = rgb(.2, .3, .3, alpha = 0.5), col = NA, cex = 1.5,
     axes=T, main = "", 
     xlab = "Climate dissimilarity",
     ylab = "Dissimilarity aboveground plane")
mod<-lmodel2(as.dist(dissimBiomeAbove)~as.dist(dissim_biome_Clmate))
curve(expr= mod$regression.results$Intercept[2]+ mod$regression.results$Slope[2]*x,
      add=T, col=rgb(.2, .3, .3), lwd=5, lty=1)


vegan::mantel(dissimBiomeBelow, dissim_biome_Clmate)
plot(as.dist(dissim_biome_Clmate), as.dist(dissimBiomeBelow), pch = 21, 
     bg = rgb(.2, .3, .3, alpha = 0.5), col = NA, cex = 1.5,
     axes=T, main = "", 
     xlab = "Climate dissimilarity",
     ylab = "Dissimilarity fine-root plane")
mod<-lmodel2(as.dist(dissimBiomeBelow)~as.dist(dissim_biome_Clmate))
curve(expr= mod$regression.results$Intercept[2]+ mod$regression.results$Slope[2]*x,
      add=T, col=rgb(.2, .3, .3), lwd=5, lty=1)

