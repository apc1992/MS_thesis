#!/usr/bin/env Rscript

# this script includes the steps for an enrichment analysis(EA) with R packages. Here we put as example our dataset analyzed with a customized ontology based on MapMan system #

# install and load required packages #
install.packages(c('gprofiler2','plotly','clusterProfiler','string','htmlwidgets'))
x = c("gprofiler2","plotly","clusterProfiler","string","htmlwidgets")
lapply(x,require, character.only=TRUE)

### Here there is an example of a list of common upregulated genes from the different time points, other EA can be performed providing the list of genes of interest ###

# path to the list of genes #
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (list of genes is missing).", call.=FALSE)
}

path = args[1]
base=path
setwd(base)

# List of genes #
query0 = read.table("./common_upregulated.txt", header = F, sep = '\t')
query0 = query0$V1
# providing custom annotation in GMT format #
custom_id=upload_GMT_file(gmtfile = "./MapMan_vitis_VCost_LO1.gmt")
# Analysis #
gostres0 = gost(query = query0, organism = custom_id, significant = T,
               correction_method = c("fdr"), evcodes = TRUE)
# Modifying title label #
gostres0[1]$result$query = "Common Up-regulated genes"
# html plot #
P0 <- gostplot(gostres0, interactive = TRUE)
# csv result #
PT0 = gostplot(gostres0, interactive = FALSE)
# Text file generation with the results #
df_PT0 <- data.frame(PT0[1]$data)
df_PT0 <- apply(df_PT0,2,as.character)
write.csv(df_PT0, "enrich_table_upregulated.csv", row.names=FALSE)
# plotting in html format #
htmlwidgets::saveWidget(as_widget(P0), "Common Up-regulated genes.html", selfcontained = TRUE)
# plotting with highlights, bear in mind that the terms must be supplied manually from visual curation #
terminos0_list = PT0$data$term_id
terminos0 = c( "2.1.1","3", "3.1","3.1.4","4","4.1",
               "4.1.5", "9", "9.2","9.2.2","9.2.2.1.1.1",
               "9.2.2.1.1","9.2.2.1","9.2.1","9.2.1.1",
               "21","24","24.2.4","24.2.4.1","24.2.4.1.1","24.2")
# Apply visual modifications to the plot #
PT0_publish = publish_gostplot(PT0, filename = "common_upregulated_genes.pdf", width = 20, height = 15, highlight_terms = terminos0)
PT0_publish + theme(plot.title  = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"))
