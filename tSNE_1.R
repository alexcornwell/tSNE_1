
install("DepecheR")             
install.packages("Rtsne")

BiocManager::install("diffcyt")
BiocManager::install("flowCore")

library(Rtsne)
library(DepecheR)  
library(diffcyt)  
library(flowCore)  #contains functions to transform data

#import subsetted sample data 
data_cytof <- read.csv(file = "./sample_subset.csv")

#remove duplicate rows from data, and remove unnecessary columns 
data_unique <- unique(data_cytof[, 3:63])  

#mass cytometry data should be transformed prior to clustering. 
#The raw data follows an approximately log-normal distribution. Transforming 
#with arcsinh improves clustering performance and allows positive and 
#negative populations to be distinguished more clearly. 
#(https://rdrr.io/bioc/diffcyt/man/transformData.html)

#define arcsinh function to transform data
#A cofactor of 5 is recomended for cyTOF data (Bendall et al. (2011), Science)
asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=1, b=5, c=1)

#apply transformation
data_cytof2 <- asinhTrans(data_unique)

#use this package to get a neat summary heatmap of parameters
cytofD <- depeche(data_cytof2)

#this function generates a tSNE, the perplexity can be changed to control,
#local vs global relationships, (e.g. lower values may bias local relationships
#when mapping clusters)
cytofSNE <- Rtsne(data_cytof2, pca=FALSE, 
                    check_duplicates = FALSE, perplexity = 30, max_iter = 2000)
  
#tSNE plot of the whole population using Base R plot
plot(cytofSNE$Y,col=cytofD$clusterVector, asp=1, cex = 0.05, xlab = "", 
       ylab = "")

#better looking plot package for tSNEs that automatically exports, high res
#.png file
dColorPlot(colorData = cytofD$clusterVector, xYData = cytofSNE$Y, 
           colorScale = "viridis", densContour = FALSE, plotName = "All_data")

#I defined these populations on web @ floreada.io with RAW concatenated .fcs files.
#The two days (8 & 11) are from the same treatment group. 
#I concatenated all events prior to phenotyping to minimize any bias due to 
#difference in expression levels between the days. The populations were gated for specific 
#phenotype and the events were exported as .csv and concatenated in BASH script with 
#populations named in "relative.event" (column 2).The Absolute.event column one 
#indicating the day (column 1)

#now I will highlight the events to determine if a phenotype corresponds to
#a cluster

plot(cytofSNE$Y[data_cytof$Relative.Event=="Bcell",], col=cytofD$clusterVector, asp=1, 
                cex = 0.05, xlab = "", ylab = "")
plot(cytofSNE$Y[data_cytof$Relative.Event=="CTL",], col=cytofD$clusterVector, asp=1, 
     cex = 0.05, xlab = "", ylab = "")
plot(cytofSNE$Y[data_cytof$Relative.Event=="Treg",], col=cytofD$clusterVector, asp=1, 
     cex = 0.05, xlab = "", ylab = "")
plot(cytofSNE$Y[data_cytof$Relative.Event=="TH",], col=cytofD$clusterVector, asp=1, 
     cex = 0.05, xlab = "", ylab = "")

#It is apparant that Tcells and Bcells are identified as seperate clusters 
#within the tSNE. tSNE is interesting in that these populations were grouped 
#algorhythmically and relationships between them were formed without 
#subjective bias.

#how do the two day 8 vs 11 compare?
plot(cytofSNE$Y[data_cytof$Absolute.Event==8,], col=cytofD$clusterVector, asp=1, 
     cex = 0.05, xlab = "", ylab = "")
plot(cytofSNE$Y[data_cytof$Absolute.Event==11,], col=cytofD$clusterVector, asp=1, 
     cex = 0.05, xlab = "", ylab = "")

#It is apparant that Day 8 data correpond to elevated Tcell activities especially
#cytotoxic t cells are very much higher. At,
#Day 11 there appears to be elevated adaptive immune stimulation, and increased
#plasma cell (Mature Bcell) formation

#now use this package to create nice quality .png plots

#define the density contour 
densContour <- dContours(cytofSNE$Y)

dDensityPlot(xYData = cytofSNE$Y[data_cytof$Relative.Event=="TH",], 
             plotName = "Thelper", colorScale="red", densContour = densContour)
dDensityPlot(xYData = cytofSNE$Y[data_cytof$Relative.Event=="Treg",], 
             plotName = "Treg", colorScale="orange", densContour = densContour)
dDensityPlot(xYData = cytofSNE$Y[data_cytof$Relative.Event=="CTL",], 
             plotName = "CTL", colorScale="yellow", densContour = densContour)
dDensityPlot(xYData = cytofSNE$Y[data_cytof$Relative.Event=="Bcell",], 
             plotName = "Bcell", colorScale="green", densContour = densContour)
dDensityPlot(xYData = cytofSNE$Y[data_cytof$Relative.Event>=1,], 
             plotName = "Lymphocyte", colorScale="blue", densContour = densContour)
dDensityPlot(xYData = cytofSNE$Y[data_cytof$Absolute.Event==8,], 
             plotName = "Day8", colorScale="green", densContour = densContour)
dDensityPlot(xYData = cytofSNE$Y[data_cytof$Absolute.Event==11,], 
             plotName = "Day11", colorScale="blue", densContour = densContour)

#doesnt work so well bc cytof data negative-peak is so dominant. 
dViolins(cytofD$clusterVector, inDataFrame = data_cytof2, 
         plotClusters = 1, plotElements = cytofD$essenceElementList)

sessionInfo()
