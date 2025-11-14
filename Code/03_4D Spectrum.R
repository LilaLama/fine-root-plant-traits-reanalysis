### In this script we plot the functional space based on species with fully empirical information for both aboveground AND fine-root traits (301 species). The resulting space is four dimensional, so that we represent plots for each pair of dimensions (Extended Data). Note that Fig. 1 in the main text represents a subset of this, showing only the aboveground (C1-C2) and fine-root (C3-C4) planes.
cat(paste0("\n\n Starting script #3 \n\n"))

# Export to Figures
png(filename = "Figures/Spectra_Above_Below_combined_full.png",
        width = 10, height = 12, units = "in", res = 300)

colGradient <- c("white",  "yellow", "red")
gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
ncolors <- 1000
ColorRamp <- rev(gradientColorsF(ncolors))
contourLevels <- c(0.5, 0.99)
cexMain <- 1.25
thickCountour <- c(4, 3)

data <- readRDS("data/PCATotal_CompleteObs.rds")

layout(matrix(c(1, 0, 0, 7, 
                2, 4, 0, 7,
                3, 5, 6, 7), nrow=3, byrow=T), 
       widths = c(0.3, 0.3, 0.3, 0.1))
###### FULL SPECTRUM (ALL TRAITS), WITH THREE DIMENSIONS:
### ALL SPECTRA FIRST TWO AXES:
par(mar = c(5, 5, 2, 0.3), 
    oma = c(1, 1, 1, 1), 
    mgp = c(2.2, 0.3, 0), 
    cex.axis = 1.5, cex.lab = 1.8)
dimensions <- data$dimensions
TPDsUse <- data$TPDs2D
for(i in 1:length(TPDsUse)){
  TPDsAux <- TPDsUse[[i]]
  traitsUSE <- data$traitsUse2D[[i]]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(traitsUSE),
                       dimnames = list(commNames, rownames(traitsUSE)))
  commMatrix["All", ] <- 1
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab <- vector("expression", 1)
  if(i > 5)   xlab[1]<- substitute(expression(paste("C3"["F"]," (", MYVALUE, "%)")), 
                                   list( MYVALUE = round(100 * data$Variance[3], 2)))[2]
  if(i <= 5)   xlab[1]<- substitute(expression(paste("C2"["F"]," (", MYVALUE, "%)")), 
                                    list( MYVALUE = round(100 * data$Variance[2], 2)))[2]
  if(i <= 3)   xlab[1]<- substitute(expression(paste("C1"["F"]," (", MYVALUE, "%)")), 
                                    list( MYVALUE = round(100 * data$Variance[1], 2)))[2]
  
  ylab<-vector("expression",1)
  ylab[1]<- substitute(expression(paste("C4"["F"]," (", MYVALUE, "%)")), 
                                   list( MYVALUE = round(100 * data$Variance[4], 2)))[2]
  if(i == 1)   ylab[1]<- substitute(expression(paste("C2"["F"]," (", MYVALUE, "%)")), 
                                    list( MYVALUE = round(100 * data$Variance[2], 2)))[2]
  if(i == 2 | i == 4)   ylab[1]<- substitute(expression(paste("C3"["F"]," (", MYVALUE, "%)")), 
                                    list( MYVALUE = round(100 * data$Variance[3], 2)))[2]
  
  imageMat <- imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[, 1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[, 2])
  limX <- c(-5, 5)
  limY <- c(-5, 5)
  image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1)
  points(traitsUSE, pch=21, cex=0.2, bg = rgb(.5, .5, .5, alpha = 0.5), col=NA)
  box(which="plot")
  axis(1,tcl=0.3,lwd=0.8)
  axis(2, las=1, tcl=0.3,lwd=0.8)
  contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.1)),
          drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey40", add=T)
  contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=contourLevels,
          drawlabels = T, labcex = 0.8, lwd=thickCountour, lty=1, col="black", add=T)
}

###GRADIENT LEGEND:
par(mar=c(1,1.5,2,1))
legend_image <- as.raster(matrix(rev(ColorRamp), ncol=1))
plot(c(0, 2),c(-0.1, 1),type = 'n', axes = F,xlab = '', ylab = '')
mtext(text = 'Species density\nquantiles', side=3,  cex = 0.9 * cexMain, line=-1)
text(x=1.5, y =c(0, 0.5, 0.99), labels = c(0, 0.5, 0.99), cex=0.9 * cexMain)
rasterImage(legend_image, xleft=0, ybottom=0, xright=1, ytop=1)
heights <- c(0, 0.5, 0.99)
for(i in 1:length(heights)){
  lines(x=c(0,1), y=rep(heights[i],2), lwd=0.5)
}

dev.off()
