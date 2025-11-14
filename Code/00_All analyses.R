### Executing the following code in order performs the main analyses and figures presented in the manuscript "Vascular plant trait diversity is globally higher above- than belowground" by Carmona et al. Before running it, the working directory must be set so that it contains the different .R files, and the data files provided are within a subfolder called "data"
# After each script we show the execution time in a Windows 10 PC (Intel i7-8700 CPU 3.2GHz, RAM 32GB. R version 3.6.0), for reference
# note that a series of packages have to be installed and loaded for the different scripts to run properly:
library(paran)
library(psych)
library(TPD)
library(shape) #Only used for plotting
library(RColorBrewer) #Only used for plotting
library(ade4)
library(mvtnorm)
library(missForest)
library(V.PhyloMaker)
library(vegan)
library(FD)
library(lmodel2)
library(plotbiomes)
library(plotrix)
library(here)

dir.create(file.path("Figures"), showWarnings = FALSE)
#0. Loading a set of functions to perform different analyses (versions of functions from the TPD package that allow for faster calculations, better plots, and summary of results from null models)
source("code/Aux_Functions.R")

#1. Estimations of functional spaces with complete empirical information
system.time(source("code/01_Full empirical information spectrum.R")) # 4.06 min

#1b. Estimations of functional space based on correlations (eigenanalysis)
system.time(source("code/01b_Eigenanalysis from correlation matrix.R")) # 0.13 secs

#2. Plot of the functional spaces based on species with complete empirical information aboveground (2630 species) and fine-root traits (748 species)
system.time(source("code/02_Spectra Only above-Only below.R")) # 10.6 sec

#3. Plot of the functional space based on species with complete empirical information for both aboveground  and fine-root traits (301 species)
system.time(source("code/03_4D Spectrum.R")) #6.7 sec

#4. Procrustes analyses to compare individual spaces between them and with the corresponding planes of the total space based on all 10 traits
system.time(source("code/04_Procrustes_Analyses.R")) # 5.69 sec

#5. Comparison of observed distribution of species with null models based on multivariate normal distributions
system.time(source("code/05_Comparing_with_multivariate_Normal.R")) # 22.8 hours

#6. Imputation of trait values for species with incomplete information and creation of functional space based on imputed dataset
system.time(source("code/06_Imputed information spectrum.R")) # 13.2 min

#7. PERMANOVA analyses
system.time(source("code/07_Permanova_Analyses.R")) # 5.6 min

#8. Dissimilarity analyses, including Fig. 2
system.time(source("code/08_Dissimilarity_Analyses.R")) # 33.5 sec

#9. Functional richness/Redundancy analyses, including Fig. 3
system.time(source("code/09_Redundancy_Analyses.R")) # 44.7 min

