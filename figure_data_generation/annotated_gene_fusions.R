rm(list=ls())
## this script is wittern to display the gene fusion statistical and annotation results
library(ggplot2)
## load the fusion results
setwd("C:/Users/ZLiu/Desktop/genefusion/Updated_fusion_results_2018.06.05")
load("Chimerascan_results.RData")
load("SOAPfuse_results.RData")
load("TopHat_results.RData")
#
frequency_fusions_Chimerascan=as.data.frame(table(Chimerascan_results[,2]))
ratio_chimerascan = length(which(frequency_fusions_Chimerascan$Freq==1))/length(frequency_fusions_Chimerascan$Freq)

frequency_fusions_SOAPfuse=as.data.frame(table(SOAPfuse_results[,2]))
ratio_SOAPfuse = length(which(frequency_fusions_SOAPfuse$Freq==1))/length(frequency_fusions_SOAPfuse$Freq)

frequency_fusions_TopHat=as.data.frame(table(TopHat_results[,2]))
ratio_TopHat = length(which(frequency_fusions_TopHat$Freq==1))/length(frequency_fusions_TopHat$Freq)
#
frequency_report_fusions_Chimerascan = as.data.frame(table(Chimerascan_results[which(Chimerascan_results[,9]==1),2]))
kk = sort(frequency_report_fusions_Chimerascan$Freq,decreasing = TRUE,index.return = TRUE)
ranked_report_Chimerascan = frequency_report_fusions_Chimerascan$Var1[kk$ix]

gene_list_report_Chimerascan = union(Chimerascan_results[which(Chimerascan_results[,9]==1),3],Chimerascan_results[which(Chimerascan_results[,9]==1),5])
write.table(gene_list_report_Chimerascan,file='sky1.txt',col.names=F,row.names = F)

frequency_report_fusions_SOAPfuse = as.data.frame(table(SOAPfuse_results[which(SOAPfuse_results[,9]==1),2]))
kk = sort(frequency_report_fusions_SOAPfuse$Freq,decreasing = TRUE,index.return = TRUE)
ranked_report_SOAPfuse = frequency_report_fusions_SOAPfuse$Var1[kk$ix]

gene_list_report_SOAPfuse = union(SOAPfuse_results[which(SOAPfuse_results[,9]==1),3],SOAPfuse_results[which(SOAPfuse_results[,9]==1),5])
write.table(gene_list_report_SOAPfuse,file='sky2.txt',col.names=F,row.names = F)

frequency_report_fusions_TopHat = as.data.frame(table(TopHat_results[which(TopHat_results[,9]==1),2]))
kk = sort(frequency_report_fusions_TopHat$Freq,decreasing = TRUE,index.return = TRUE)
ranked_report_TopHat = frequency_report_fusions_TopHat$Var1[kk$ix]

gene_list_report_TopHat = union(TopHat_results[which(TopHat_results[,9]==1),3],TopHat_results[which(TopHat_results[,9]==1),5])
write.table(gene_list_report_TopHat,file='sky3.txt',col.names=F,row.names = F)
#
