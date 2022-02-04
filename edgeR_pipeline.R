#!/usr/bin/env Rscript

# this script includes the differential expression analysis(DEA) with edgeR package and the results displayment for the 4 time points that conform our experimental design #

# install packages #
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HTSFilter")
BiocManager::install("edgeR")
install.packages(c('ggplot2','cluster','ggfortify','ggplotly','plotly','ggpubr','ggrepel'))

# load packages #
x = c("HTSFilter","edgeR","ggplot2","cluster","ggfortify","plotly","ggpubr","ggrepel")
lapply(x,require, character.only=TRUE)

# path to the row count matrix and the metadata table #
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (path to the folder that contains the raw counts and metadata).", call.=FALSE)
}

path = args[1]
base=path
setwd(base)

# import counts matrix and metadata dataframe #
counts = read.table("./corrected_counts_3.25.txt",sep="\t",header=T,row.names=1)
metadata = read.table("./metadata.txt",sep="\t",header=T)

### NORMALIZATION ###

# merge experimental conditions in an additional metadata column, factor #
metadata$Factor = as.factor(paste(metadata$Colour,"_",metadata$Time,sep = ""))
# apply the filtering + TMM normalization #
filtro = HTSFilter(counts,metadata$Factor,normalization="TMM")
# extract the filtered counts into a new dataframe #
counts_filtrados = as.data.frame(filtro$filteredData)

## PCA plot ### 

pca_all_samples = prcomp(t(counts_filtrados))
graf = autoplot(fanny(pca_all_samples$x, 4), frame = TRUE,label= TRUE, shape= FALSE, frame.type = 'norm')
# save plot in PDF/PNG format #
ggsave("PCA_global.png",device = "png",width = 8,height = 7,dpi = 600, path = base)
ggsave("PCA_global.pdf",device = "pdf",width = 8,height = 7,dpi = 600, path = base)
# create and save the interactive version of the plot in html #
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"pca.html")
# customize the original plot #
g = graf + theme_bw()
gg = g + theme(axis.title.x = element_text(size = 18, face = "bold"), 
                 axis.title.y = element_text(size = 18, face = "bold"),
                 legend.position = "none")
gg 


### GLOBAL DEA WITH EDGER ###

# calculate the DGE list and define the model depending on experimental variables #
analysis = DGEList(counts = counts_filtrados)
modelo = with(metadata,model.matrix(~ 0 + Colour))
colnames(modelo) = c("red","white")
# apply contrasts depending on the defined model #
contrast = makeContrasts(White_vs_Red = white - red,levels = modelo)
analysis = calcNormFactors(analysis,method = "TMM")
analysis = estimateDisp(analysis,design = modelo)
fit = glmFit(analysis,design = modelo)
lrt = glmLRT(fit,contrast = contrast[,1])
# extract results #
resultados = topTags(lrt,n=Inf)[[1]]
# filter only significant results #
resultados_sig = resultados[which(resultados$FDR<=0.05),]

### T4 DEA WITH EDGER ###

# filter data corresponding to t4 #
counts_t4 = counts_filtrados[,which(metadata$Time==4)]
metadata_t4 = metadata[which(metadata$Time==4),]

# apply DEA only for t4 #
analysis = DGEList(counts = counts_t4 )
model_t4 = with(metadata_t4,model.matrix(~ 0 + Colour))
colnames(model_t4) = c("red","white")
contrast_t4 = makeContrasts(White_vs_Red = white - red,levels = model_t4)
analysis = calcNormFactors(analysis,method = "TMM")
analysis = estimateDisp(analysis,design = model_t4)
fit_t4 = glmFit(analysis,design = model_t4)
lrt_t4 = glmLRT(fit_t4,contrast = contrast_t4[,1])
resultados_t4 = topTags(lrt_t4,n=Inf)[[1]]
resultados_t4_sig = resultados_t4[which(resultados_t4$FDR<=0.05),]
# split significant outcomes into up/downregulated and store the results #
upregulated_t4 = resultados_t4[which(resultados_t4$logFC> 0.5 & resultados_t4$FDR<0.05),]
downregulated_t4 = resultados_t4[which(resultados_t4$logFC< -0.5 & resultados_t4$FDR<0.05),]
write.table(x = upregulated_t4, file = "./9WAV_up.txt", sep = "\t")
write.table(x = downregulated_t4, file = "./9WAV_down.txt", sep = "\t")

## Volcano-plot ##

# add status column in the results dataframe and update it depending on the state #
resultados_t4$status = "Not significant"
resultados_t4$status[resultados_t4$logFC> 0.5 & resultados_t4$FDR<0.05] = "Up-regulated"
resultados_t4$status[resultados_t4$logFC< -0.5 & resultados_t4$FDR<0.05] = "Down-regulated"
# generate the plot #
g=ggplot(data=resultados_t4,aes(x=logFC,y=-log10(FDR), col=status)) + 
  geom_point() + xlim(-8,8) + theme_minimal()
volcano_t4 = g + geom_vline(xintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") + 
  geom_hline(yintercept = -log10(0.05),col="blue", linetype="dotted")
# save in pdf/png format #
ggsave("volcano.png",device = "png",width = 8,height = 7,dpi = 600, path = base)
ggsave("volcano.pdf",device = "pdf",width = 8,height = 7,dpi = 600, path = base)
# generate and store the interactive version in html format #
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"volcano.html")

## MA-plot ##

# same procedure as in volcano plot #
g=ggplot(data=resultados_t4,aes(x=logCPM,y=logFC, col=status)) +
  geom_point()  + theme_minimal() + ylim(-8,8)
ma_t4 = g + geom_hline(yintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") 
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"MA.html")

  
## bar-plot ##

# generate the dataframe that contains the data of interest for the plot #
df = as.data.frame(nrow(counts_t4))
df$sig=nrow(resultados_t4_sig)
df$up= length(which(resultados_t4_sig$logFC>0.5))
df$down= length(which(resultados_t4_sig$logFC<0.5))
colnames(df)= c("Total","Significant","Up-regulated","Down-regulated")
x=colnames(df)
df=as.data.frame(t(df))
conts=df$V1
df = data.frame(x=x, counts=conts)
# display the plot #
bar_t4 = ggplot(data=df,aes(x=reorder(x, counts),y=counts)) + geom_bar(stat="identity", fill="steelblue", width = 0.5) +
  theme_minimal() + ylab("Number of Genes") + xlab("") + coord_flip()
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"barplot.html")


## THE CODE BELLOW HAS BEEN APPLIED TO COMBINE THE 3 PLOTS CORRESPONDING TO 9WAV(T4), THE SAME PROCEDURE HAS BEEN APPLIED TO ANY COMBINATION OF PLOTS  ##

t4_combined = ggarrange(bar_t4, ggarrange(volcano_t4,ma_t4, labels= c("Volcano","MA"), common.legend = T, 
                        legend = "right"), nrow=2, labels= "Summary")
t4_combined

### T3 DEA WITH EDGER ###

# The same protocol as in t4 has been applied to all the time points #

counts_t3 = counts_filtrados[,which(metadata$Time==3)]
metadata_t3 = metadata[which(metadata$Time==3),]

# DEA #
analysis = DGEList(counts = counts_t3 )
model_t3 = with(metadata_t3,model.matrix(~ 0 + Colour))
colnames(model_t3) = c("red","white")
contrast_t3 = makeContrasts(White_vs_Red = white - red,levels = model_t3)
analysis = calcNormFactors(analysis,method = "TMM")
analysis = estimateDisp(analysis,design = model_t3)
fit_t3 = glmFit(analysis,design = model_t3)
lrt_t3 = glmLRT(fit_t3,contrast = contrast_t3[,1])
resultados_t3 = topTags(lrt_t3,n=Inf)[[1]]
resultados_t3_sig = resultados_t3[which(resultados_t3$FDR<=0.05),]
upregulated_t3 = resultados_t3[which(resultados_t3$logFC> 0.5 & resultados_t3$FDR<0.05),]
downregulated_t3 = resultados_t3[which(resultados_t3$logFC< -0.5 & resultados_t3$FDR<0.05),]
write.table(x = upregulated_t3, file = "./6WAV_up.txt", sep = "\t")
write.table(x = downregulated_t3, file = "./6WAV_down.txt", sep = "\t")

##Volcano-plot##

resultados_t3$status = "Not significant"
resultados_t3$status[resultados_t3$logFC> 0.5 & resultados_t3$FDR<0.05] = "Up-regulated"
resultados_t3$status[resultados_t3$logFC< -0.5 & resultados_t3$FDR<0.05] = "Down-regulated"
g=ggplot(data=resultados_t3,aes(x=logFC,y=-log10(FDR), col=status)) + 
  geom_point() + xlim(-8,8) + theme_minimal()
volcano_t3 = g + geom_vline(xintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") + 
  geom_hline(yintercept = -log10(0.05),col="blue", linetype="dotted")
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"volcano.html")


##MA-plot##

g=ggplot(data=resultados_t3,aes(x=logCPM,y=logFC, col=status)) +
  geom_point()  + theme_minimal() + ylim(-8,8)
ma_t3 = g + geom_hline(yintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") 
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"MA.html")

##bar-plot##
df = as.data.frame(nrow(counts_t3))
df$sig=nrow(resultados_t3_sig)
df$up= length(which(resultados_t3_sig$logFC>0.5))
df$down= length(which(resultados_t3_sig$logFC<0.5))
colnames(df)= c("Total","Significant","Up-regulated","Down-regulated")
x=colnames(df)
df=as.data.frame(t(df))
conts=df$V1
df = data.frame(x=x, counts=conts)

bar_t3= ggplot(data=df,aes(x=reorder(x, counts),y=counts)) + geom_bar(stat="identity", fill="steelblue", width = 0.5) +
  theme_minimal() + ylab("Number of Genes") + xlab("") + coord_flip()
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"barplot.html")

### T2 DEA WITH EDGER ###

counts_t2 = counts_filtrados[,which(metadata$Time==2)]
metadata_t2 = metadata[which(metadata$Time==2),]

# DEA #
analysis = DGEList(counts = counts_t2 )
model_t2 = with(metadata_t2,model.matrix(~ 0 + Colour))
colnames(model_t2) = c("red","white")
contrast_t2 = makeContrasts(White_vs_Red = white - red,levels = model_t2)
analysis = calcNormFactors(analysis,method = "TMM")
analysis = estimateDisp(analysis,design = model_t2)
fit_t2 = glmFit(analysis,design = model_t2)
lrt_t2 = glmLRT(fit_t2,contrast = contrast_t2[,1])
resultados_t2 = topTags(lrt_t2,n=Inf)[[1]]
resultados_t2_sig = resultados_t2[which(resultados_t2$FDR<=0.05),]
upregulated_t2 = resultados_t2[which(resultados_t2$logFC> 0.5 & resultados_t2$FDR<0.05),]
downregulated_t2 = resultados_t2[which(resultados_t2$logFC< -0.5 & resultados_t2$FDR<0.05),]
write.table(x = upregulated_t2, file = "./3WAV_up.txt", sep = "\t")
write.table(x = downregulated_t2, file = "./3WAV_down.txt", sep = "\t")

## Volcanot-plot ##

resultados_t2$status = "Not significant"
resultados_t2$status[resultados_t2$logFC> 0.5 & resultados_t2$FDR<0.05] = "Up-regulated"
resultados_t2$status[resultados_t2$logFC< -0.5 & resultados_t2$FDR<0.05] = "Down-regulated"
g=ggplot(data=resultados_t2,aes(x=logFC,y=-log10(FDR), col=status)) + 
  geom_point() + xlim(-8,8) + theme_minimal()
volcano_t2 = g + geom_vline(xintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") + 
  geom_hline(yintercept = -log10(0.05),col="blue", linetype="dotted")

volcano_t2 = volcano_t2 + theme(axis.title.x = element_text(size = 18, face = "bold"), 
                   axis.title.y = element_text(size = 18, face = "bold"),
                    legend.title = element_blank())
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"volcano.html")

##MA-plot##

g=ggplot(data=resultados_t2,aes(x=logCPM,y=logFC, col=status)) +
  geom_point()  + theme_minimal() + ylim(-8,8)
ma_t2 = g + geom_hline(yintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") 
ma_t2 = ma_t2 + theme(axis.title.x = element_text(size = 18, face = "bold"), 
              axis.title.y = element_text(size = 18, face = "bold"),
              legend.title = element_blank())
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"MA.html")

##bar-plot##

df = as.data.frame(nrow(counts_t2))
df$sig=nrow(resultados_t2_sig)
df$up= length(which(resultados_t2_sig$logFC>0.5))
df$down= length(which(resultados_t2_sig$logFC<0.5))
colnames(df)= c("Total","Significant","Up-regulated","Down-regulated")
x=colnames(df)
df=as.data.frame(t(df))
conts=df$V1
df = data.frame(x=x, counts=conts)

bar_t2 = ggplot(data=df,aes(x=reorder(x, counts),y=counts)) + geom_bar(stat="identity", fill="steelblue", width = 0.5) +
  theme_minimal() + ylab("Number of Genes") + xlab("") + coord_flip()
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"barplot.html")  

### T1 DEA WITH EDGER ###

counts_t1 = counts_filtrados[,which(metadata$Time==1)]
metadata_t1 = metadata[which(metadata$Time==1),]

# DEA #
analysis = DGEList(counts = counts_t1 )
model_t1 = with(metadata_t1,model.matrix(~ 0 + Colour))
colnames(model_t1) = c("red","white")
contrast_t1 = makeContrasts(White_vs_Red = white - red,levels = model_t1)
analysis = calcNormFactors(analysis,method = "TMM")
analysis = estimateDisp(analysis,design = model_t1)
fit_t1 = glmFit(analysis,design = model_t1)
lrt_t1 = glmLRT(fit_t1,contrast = contrast_t1[,1])
resultados_t1 = topTags(lrt_t1,n=Inf)[[1]]
resultados_t1_sig = resultados_t1[which(resultados_t1$FDR<=0.05),]
upregulated_t1 = resultados_t1[which(resultados_t1$logFC> 0.5 & resultados_t1$FDR<0.05),]
downregulated_t1 = resultados_t1[which(resultados_t1$logFC< -0.5 & resultados_t1$FDR<0.05),]
write.table(x = upregulated_t1,file="./0WAV_up.txt", sep = "\t")
write.table(x = downregulated_t1,file = "./0WAV_down.txt", sep = "\t")

##Volcano-plot##

resultados_t1$status = "Not significant"
resultados_t1$status[resultados_t1$logFC> 0.5 & resultados_t1$FDR<0.05] = "Up-regulated"
resultados_t1$status[resultados_t1$logFC< -0.5 & resultados_t1$FDR<0.05] = "Down-regulated"
g=ggplot(data=resultados_t1,aes(x=logFC,y=-log10(FDR), col=status)) + 
  geom_point() + xlim(-8,8) + theme_minimal()
volcano_t1 = g + geom_vline(xintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") + 
  geom_hline(yintercept = -log10(0.05),col="blue", linetype="dotted")
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"volcano.html")

##MA-plot##

g=ggplot(data=resultados_t1,aes(x=logCPM,y=logFC, col=status)) +
  geom_point()  + theme_minimal() + ylim(-8,8)
ma_t1 = g + geom_hline(yintercept = c(-0.5, 0.5),col = "blue", linetype="dotted") 
p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"MA.html")

##bar-plot##
df = as.data.frame(nrow(counts_t1))
df$sig=nrow(resultados_t1_sig)
df$up= length(which(resultados_t1_sig$logFC>0.5))
df$down= length(which(resultados_t1_sig$logFC<0.5))
colnames(df)= c("Total","Significant","Up-regulated","Down-regulated")
x=colnames(df)
df=as.data.frame(t(df))
conts=df$V1
df = data.frame(x=x, counts=conts)

bar_t1 = ggplot(data=df,aes(x=reorder(x, counts),y=counts)) + geom_bar(stat="identity", fill="steelblue", width = 0.5) +
  theme_minimal() + ylab("Number of Genes") + xlab("") + coord_flip()

p = ggplotly()
htmlwidgets::saveWidget(as_widget(p),"barplot.html") 
