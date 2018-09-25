rm(list=ls())
## This script is used to map the DEG from different fusion detection tools onto immune system cell lines
setwd('C:/Users/ZLiu/Desktop/genefusion/Updated_results_2018.06.13')
##
ChimeraScan_DEG = read.table('ChimeraScan_DEG.csv',header=T,sep=',')
SOAPfuse_DEG = read.table('SOAPfuse_DEG.csv',header=T,sep=',')
TopHat_DEG = read.table('TopHat_DEG.csv',header=T,sep=',')

## map gene symbols to entrz Gene ID
library("org.Hs.eg.db")
up_genes=as.character(ChimeraScan_DEG[c(1:500),1])
down_genes = as.character(ChimeraScan_DEG[c(501:1000),1])
genes_chimeraScan_up = as.numeric(unlist(mget(x=up_genes,envir=org.Hs.egALIAS2EG)))
genes_chimeraScan_down = as.numeric(unlist(mget(x=down_genes,envir=org.Hs.egALIAS2EG)))
rm(up_genes,down_genes)
#
up_genes=as.character(SOAPfuse_DEG[c(1:500),1])
down_genes = as.character(SOAPfuse_DEG[c(501:1000),1])
genes_SOAPfuse_up = as.numeric(unlist(mget(x=up_genes,envir=org.Hs.egALIAS2EG)))
genes_SOAPfuse_down = as.numeric(unlist(mget(x=down_genes,envir=org.Hs.egALIAS2EG)))
rm(up_genes,down_genes)
#
up_genes=as.character(TopHat_DEG[c(1:500),1])
down_genes = as.character(TopHat_DEG[c(501:1000),1])
genes_TopHat_up = unlist(mget(x=up_genes,envir=org.Hs.egALIAS2EG))
genes_TopHat_down = unlist(mget(x=down_genes,envir=org.Hs.egALIAS2EG))
rm(up_genes,down_genes)

## load immune cell line signatures from ImmGen
immune_signatures = read.table('immune_data.csv',header=T,sep=',')
genes_entrz_immune = immune_signatures[,1]
immune_signatures = immune_signatures[,-1]
## similarity betwwen DEGs
simiarlity_ChimeraScan = matrix(0,ncol(immune_signatures),1)
simiarlity_SOAPfuse = matrix(0,ncol(immune_signatures),1)
simiarlity_TopHat = matrix(0,ncol(immune_signatures),1)
for (i in 1:ncol(immune_signatures)){
  index = sort(immune_signatures[,i],decreasing = T, index.return = T)
  ranked_genes_entrz_immune = genes_entrz_immune[index$ix]
  up_immune = ranked_genes_entrz_immune[1:500]
  down_immune = tail(ranked_genes_entrz_immune,500)
  simiarlity_ChimeraScan[i]= (length(intersect(genes_chimeraScan_up,up_immune))+length(intersect(genes_chimeraScan_down,down_immune)))/length(union(c(genes_chimeraScan_up,up_immune),c(genes_chimeraScan_down,down_immune)))
  simiarlity_SOAPfuse[i]= (length(intersect(genes_SOAPfuse_up,up_immune))+length(intersect(genes_SOAPfuse_down,down_immune)))/length(union(c(genes_SOAPfuse_up,up_immune),c(genes_SOAPfuse_down,down_immune)))
  simiarlity_TopHat[i]= (length(intersect(genes_TopHat_up,up_immune))+length(intersect(genes_TopHat_down,down_immune)))/length(union(c(genes_TopHat_up,up_immune),c(genes_TopHat_down,down_immune)))
}
## ranked order the immmune cell lines
cell_line_names = colnames(head(immune_signatures,0))
ranked_similarity_ChimeraScan = sort(simiarlity_ChimeraScan,decreasing = T, index.return = T)
ranked_similarity_SOAPfuse = sort(simiarlity_SOAPfuse,decreasing = T, index.return = T)
ranked_similarity_TopHat = sort(simiarlity_TopHat,decreasing = T, index.return = T)
ranked_cellline_ChimeraScan = cell_line_names[ranked_similarity_ChimeraScan$ix]
ranked_cellline_SOAPfuse = cell_line_names[ranked_similarity_SOAPfuse$ix]
ranked_cellline_TopHat = cell_line_names[ranked_similarity_TopHat$ix]

##
POP_CS = matrix(0,nrow(simiarlity_ChimeraScan),1)
POP_CT = matrix(0,nrow(simiarlity_ChimeraScan),1)
POP_ST = matrix(0,nrow(simiarlity_ChimeraScan),1)
for (i in 1:nrow(simiarlity_ChimeraScan)){
  POP_CS[i] = length(intersect(ranked_cellline_ChimeraScan[c(1:i)],ranked_cellline_SOAPfuse[c(1:i)]))/i
  POP_CT[i] = length(intersect(ranked_cellline_ChimeraScan[c(1:i)],ranked_cellline_TopHat[c(1:i)]))/i
  POP_ST[i] = length(intersect(ranked_cellline_SOAPfuse[c(1:i)],ranked_cellline_TopHat[c(1:i)]))/i
}
pop_matrix = cbind(POP_CS,POP_CT,POP_ST)
top_ten_cell_types = cbind(ranked_cellline_ChimeraScan[c(1:304)],ranked_similarity_ChimeraScan$x,
                           ranked_cellline_SOAPfuse[c(1:304)],ranked_similarity_SOAPfuse$x,
                           ranked_cellline_TopHat[c(1:304)],ranked_similarity_TopHat$x)
write.table(top_ten_cell_types,file='top_ten_immune_celllines.txt')
write.table(pop_matrix,file='POP.txt',col.names=F)
