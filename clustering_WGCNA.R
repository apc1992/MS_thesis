#!/usr/bin/env Rscript

# this script contains the clustering algorithm applied to the time-course DEA results with WGCNA R package. Parameters used have been customized and are not valid for general usage #

# path to the raw counts of the genes established as significant with the time-course pipeline  #
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (counts matrix with length is missing).", call.=FALSE)
}

path = args[1]
base=path
setwd(base)

# input read counts must be normalized with TPM #
# from DEA results we need to extract the raw counts and the transcript length column #
# divide raw counts by the corresponding transcript length #
rpk=counts_length[,2:25]/counts_length$Length
# calculate the scaling factor #
scalingFactor=colSums(rpk)/1000000
# converting to the required dataframe format #
tpm=t(t(rpk))/scalingFactor
tpm=as.data.frame(tpm)
# extracting only genes significant in ImpulseDE2 results #
tpm_significant=tpm[which(rownames(tpm)%in%results_DE_TC_sig$Gene),]  
# writing the results in a csv file # 
write.table(tpm_significant,file="tpm_significant.csv",sep="\t",row.names = T,col.names = T)

# load and install library 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("WGCNA")
library('WGCNA')
# setting up the input dataframe #
options(stringsAsFactors = FALSE)
vitData=read.csv(file = 'tpm_significant.csv',
                 sep = '\t', dec = '.')
row.names(vitData) = vitData$Vcost
vitData$Vcost= NULL

# Constructing  a  weighted  gene  network  entails  the  choice  of  the  soft  thresholding  powerβto which  co-expressionsimilarity is raised to calculate adjacency #
# Choose a set of soft-thresholding powers, we can construct signed or unsigned networks, in our case we are interested in identifying positive and negative correlations so signed correlations will be applied #

# Plot the results of powers graph and connectivity #
powers = c(1:40)
# Call the network topology analysis function
sft = pickSoftThreshold(t(vitData), powerVector = powers, verbose = 5, 
                        networkType = 'signed')
# Plot the results #
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power #
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
     
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Assesing scale free network topology #

datExpr = t(vitData)
ADJ1=abs(cor(datExpr,use="p"))^16
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
k=softConnectivity(datE=datExpr,power=16,
                   type = "signed")
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# constructing the network with our custom parameters, power=16 #

net_p16 = blockwiseModules(numericLabels = T,datExpr = t(vitData), power = 16,
                           networkType = "signed",
                           maxBlockSize = 12000,
                           minModuleSize = 40,
                           deepSplit = 4,
                           mergeCutHeight = 0.1)

# We extracted the list of genes that belongs to each module in a text format file #
lista = net_p16$colors
write.table(x = lista, file = 'List_clusters_genes.txt', sep = '\t', quote = F)
# To know how many clusters we have constructed #
length(table(net_p16$colors))
# To display how many genes are associated to each cluster #
table(net_p16$colors)

# dendogram construction #

mergedColors = labels2colors(net_p10$colors)
plotDendroAndColors(net_p10$dendrograms[[1]], mergedColors[net_p10$blockGenes[[1]]], "Module colors", dendroLabels = F, hang = 0.03, addGuide = TRUE, guideHang = 0.05)                           
