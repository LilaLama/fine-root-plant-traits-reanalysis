### In this script we plot the aboveground and fine-root functional spaces based on species with fully empirical information for aboveground or fine-root traits. These spaces are hence analogous to the ones published by Diaz et al (2016) for aboveground traits and Bergmann et al. (2020) for fine-root traits. For each of these spaces we represent the density of species (trait probability density) 
cat(paste0("\n\n Starting script #2 \n\n"))

# Export to Figures
png(filename = "Figures/Spectra_Above_Below_full.png",
        width = 10, height = 5, units = "in", res = 300)

colGradient <- c("white",  "yellow", "red")
gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
ncolors <- 1000
ColorRamp <- rev(gradientColorsF(ncolors))
contourLevels <- c(0.5, 0.99)
cexMain <- 1.25
thickCountour <- c(4, 3)

pcaALL <- list()
pcaALL[[1]] <- readRDS("data/PCAAboveONLY_CompleteObs.rds")
pcaALL[[2]] <- readRDS("data/PCABelowONLY_CompleteObs.rds")
names(pcaALL) <- c("Above", "Below")

par(mar = c(4, 4, 2, 0.3), 
    oma = c(1, 1, 1, 1), 
    mgp = c(2.2, 0.1, 0), 
    cex.axis = 1.25, cex.lab = 1.5,
    mfrow = c(1, 2))
for(i in 1:2){
  data <- pcaALL[[i]]
  dimensions <- data$dimensions
  traitsUSE <- data$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:dimensions)
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(traitsUSE),
                       dimnames = list(commNames, rownames(traitsUSE)))
  commMatrix["All", ] <- 1
  TPDsAux <- data$TPDs
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab<-vector("expression",1)
  xlab[1]<- substitute(expression(paste("C1"["A"]," (", MYVALUE, "%)")), 
                       list( MYVALUE = round(100 * data$Variance[1], 2)))[2]
  ylab<-vector("expression",1)
  ylab[1]<- substitute(expression(paste("C2"["A"]," (", MYVALUE, "%)")), 
                       list( MYVALUE = round(100 * data$Variance[2], 2)))[2]
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  limX <- c(-5.5, 6.5)
  limY <- c(-5.5, 6.5)
  image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1)
  # points(traitsUSE, pch=21, cex=0.5, bg = rgb(.5, .5, .5, alpha = 0.5), col=NA)
  box(which="plot")
  axis(1,tcl=0.3,lwd=0.8)
  axis(2, las=1, tcl=0.3,lwd=0.8)
  contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.1)),
          drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey40", add=T)
  contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=contourLevels,
          drawlabels = T, labcex = 0.8, lwd=thickCountour, lty=1, col="black", add=T)
}

dev.off()
