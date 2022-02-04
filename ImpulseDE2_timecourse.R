#!/usr/bin/env Rscript

# this script is intended to be as a guide to use ImpulseDE2 package for time-course DEA, the filtered counts and the experimental conditions must be supplied by the user #
# experimental dataframe Conditions (line19) must be rewritten by the user #
# giving the path to the folder that contains the filtered counts #

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (path to the folder that contains filtered counts).", call.=FALSE)
}

path = args[1]
base=path
setwd(base)

# define the required dataframes #
impulse_counts = as.matrix(filtered_counts)
Conditions = data.frame(Sample=colnames(impulse_counts),Condition=c(rep("case",12),rep("control",12)),Time=c(3,rep("0",3),rep("3",2),rep("6",3),rep("9",3),rep("0",3),rep("3",3),rep("6",3),rep("9",3))

# install and load Impulse #
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ImpulseDE2")
library('ImpulseDE2')

# Run the Impulse model #
DE_TC = runImpulseDE2(matCountData = impulse_counts,dfAnnotation = Conditions,boolCaseCtrl = T, scaQThres=0.05, boolVerbose = T)

# Extracting the results from the impulse object #
results=DE_TC@dfImpulseDE2Results
# Ensuring that all the results are significant #
results_sig=results[which(results$padj<=0.05),]
# Generating the plots for the all the gene trajectories #
AllTrajectories=plotGenes(vecGeneIDs = results_sig$Gene,objectImpulseDE2 = DE_TC,boolCaseCtrl = T,boolMultiplePlotsPerPage = T,boolSimplePlot = T)

# Apply modifications to the plots and save in PDF format #
for(i in 1:length(results_sig$Gene)) {AllTrajectories[[i]] +
    xlab("Developmental stage (WAV)") + ylab("Read counts (TMM)") +
    theme(legend.background = element_rect(fill="#bdbdbd", size=0.5, linetype="solid")) +
    theme(legend.background = element_rect(fill="gray100",size=0.5, linetype="solid", colour ="gray30")) +
    scale_x_continuous(breaks = c(0,3,6,9)) +
    scale_color_manual(values = c('#008837','#bdbdbd','#7b3294'), labels = c("White skin", "Combined","Red skin")) +
    ggsave(filename = paste("myplot",i,".pdf",sep=""),plot = last_plot(),device = "pdf", path = "./trajectories/")}

