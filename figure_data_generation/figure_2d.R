rm(list=ls())
## this script is wittern to display the gene fusion statistical and annotation results
library(ggplot2)
library(VennDiagram)
## load the fusion results
setwd("C:/Users/ZLiu/Desktop/genefusion/Updated_fusion_results_2018.06.05")
#load("Chimerascan_results.RData")
#load("SOAPfuse_results.RData")
load("TopHat_results.RData")
#
#Chimerascan_results[which(TopHat_results[,10]=='2a'),10]='2'
#SOAPfuse_results[which(TopHat_results[,10]=='2b'),10]='2'
TopHat_results[which(TopHat_results[,10]=='4S'),10]='4'
#
frequency_fusions_S1=as.data.frame(table(TopHat_results[which(TopHat_results[,10]=='1'),2]))
kk = sort(frequency_fusions_S1$Freq,decreasing = TRUE,index.return = TRUE)
hF_S1=frequency_fusions_S1$Var1[kk$ix]
rm(kk)

frequency_fusions_S2=as.data.frame(table(TopHat_results[which(TopHat_results[,10]=='2'),2]))
kk = sort(frequency_fusions_S2$Freq,decreasing = TRUE,index.return = TRUE)
hF_S2=frequency_fusions_S2$Var1[kk$ix]
rm(kk)

frequency_fusions_S3=as.data.frame(table(TopHat_results[which(TopHat_results[,10]=='3'),2]))
kk = sort(frequency_fusions_S3$Freq,decreasing = TRUE,index.return = TRUE)
hF_S3=frequency_fusions_S3$Var1[kk$ix]
rm(kk)

frequency_fusions_S4=as.data.frame(table(TopHat_results[which(TopHat_results[,10]=='4'),2]))
kk = sort(frequency_fusions_S4$Freq,decreasing = TRUE,index.return = TRUE)
hF_S4=frequency_fusions_S4$Var1[kk$ix]
rm(kk)

frequency_fusions_HR=as.data.frame(table(TopHat_results[which(TopHat_results[,12]=='HR'),2]))
kk = sort(frequency_fusions_HR$Freq,decreasing = TRUE,index.return = TRUE)
hF_HR=frequency_fusions_HR$Var1[kk$ix]
rm(kk)

S1 = hF_S1[1:50]
S2 = hF_S2[1:50]
S3 = hF_S3[1:50]
S4 = hF_S4[1:50]
HR = hF_HR[1:50]

write.table(S1,file ='S1.txt',col.names=F,row.names = F)
write.table(S2,file ='S2.txt',col.names=F,row.names = F)
write.table(S3,file ='S3.txt',col.names=F,row.names = F)
write.table(S4,file ='S4.txt',col.names=F,row.names = F)
write.table(HR,file ='HR.txt',col.names=F,row.names = F)

## input the files S1 -S5 to http://bioinformatics.psb.ugent.be/webtools/Venn/ to generate venn diagrams

