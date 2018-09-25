rm(list=ls())
## this script is wittern to display the gene fusion statistical and annotation results
library(ggplot2)
library(venneuler)
## load the fusion results
setwd("C:/Users/ZLiu/Desktop/genefusion/Updated_fusion_results_2018.06.05")
load("Chimerascan_results.RData")
load("SOAPfuse_results.RData")
load("TopHat_results.RData")

## overlap between the three algorithms (Figure 2A)
overlap12 = intersect(Chimerascan_results[,2],SOAPfuse_results[,2])
overlap23 = intersect(TopHat_results[,2],SOAPfuse_results[,2])
overlap13 = intersect(Chimerascan_results[,2],TopHat_results[,2])
overlap_2 = union(union(overlap12,overlap23),overlap13)
overlap_3 = intersect(intersect(overlap12,overlap23),overlap13)

# fusion count across the stages (Figure 2B)
index1 = match(unique(TopHat_results[which(TopHat_results[,10]=='1'),2]),TopHat_results[,2])
c1=length(unique(TopHat_results[which(TopHat_results[,10]=='1'),2]))
overlap1=length(which(TopHat_results[index1,4]==TopHat_results[index1,6]))

index2 = match(unique(TopHat_results[which(TopHat_results[,10]=='2'),2]),TopHat_results[,2])
c2=length(unique(TopHat_results[which(TopHat_results[,10]=='2'),2]))
overlap2=length(which(TopHat_results[index2,4]==TopHat_results[index2,6]))

index3 = match(unique(TopHat_results[which(TopHat_results[,10]=='2a'),2]),TopHat_results[,2])
c3=length(unique(TopHat_results[which(TopHat_results[,10]=='2a'),2]))
overlap3=length(which(TopHat_results[index3,4]==TopHat_results[index3,6]))

index4 = match(unique(TopHat_results[which(TopHat_results[,10]=='2b'),2]),TopHat_results[,2])
c4=length(unique(TopHat_results[which(TopHat_results[,10]=='2b'),2]))
overlap4=length(which(TopHat_results[index4,4]==TopHat_results[index4,6]))

index5 = match(unique(TopHat_results[which(TopHat_results[,10]=='3'),2]),TopHat_results[,2])
c5=length(unique(TopHat_results[which(TopHat_results[,10]=='3'),2]))
overlap5=length(which(TopHat_results[index5,4]==TopHat_results[index5,6]))

index6 = match(unique(TopHat_results[which(TopHat_results[,10]=='4'),2]),TopHat_results[,2])
c6=length(unique(TopHat_results[which(TopHat_results[,10]=='4'),2]))
overlap6=length(which(TopHat_results[index6,4]==TopHat_results[index6,6]))

index7 = match(unique(TopHat_results[which(TopHat_results[,10]=='4S'),2]),TopHat_results[,2])
c7=length(unique(TopHat_results[which(TopHat_results[,10]=='4S'),2]))
overlap7=length(which(TopHat_results[index7,4]==TopHat_results[index7,6]))

index8 = match(unique(TopHat_results[which(TopHat_results[,12]=='HR'),2]),TopHat_results[,2])
c8=length(unique(TopHat_results[which(TopHat_results[,12]=='HR'),2]))
overlap8=length(which(TopHat_results[index8,4]==TopHat_results[index8,6]))

c2=c2+c3+c4;c3=c5;c4=c6+c7;ch=c8;
overlap2 = overlap2+overlap3+overlap4;overlap3=overlap5;overlap4=overlap6+overlap7;overlap_ch=overlap8
rm(c5,c6,c7,c8,overlap5,overlap6,overlap7,overlap8)
dis1 = c1-overlap1;dis2=c2-overlap2;dis3=c3-overlap3;dis4=c4-overlap4;dis_ch=ch-overlap_ch
rm(list = ls(pattern = "index"))

## concordance figure (Figure S1)
concordance = cbind(tani,as.character(clinical_info[,3]))
temp=tani[which(clinical_info[,5]=="HR")]
tani=c(tani,temp)
hh=as.character(rep.int('HR',length(temp)))
tag=c(as.character(clinical_info[,3]),hh)
concordance=data.frame(tani,tag)
colnames(concordance)=c('tani','Patient_classfication')
concordance$Patient_classfication=as.factor(concordance$Patient_classfication)
#concordance=concordance[-which(concordance$tani==0),]
concordance=concordance[-which(concordance$Patient_classfication=='multilok. (2)'),]
concordance[which(concordance$Patient_classfication=='2a'),2]='2';
concordance[which(concordance$Patient_classfication=='2b'),2]='2';
concordance[which(concordance$Patient_classfication=='4S'),2]='4';
p<-ggplot(concordance, aes(x=Patient_classfication, y=tani, color=Patient_classfication),geom = "ribbon") + geom_violin(trim=FALSE)+
  geom_dotplot(stackdir='down', dotsize=0.02) + geom_jitter(shape=16, position=position_jitter(0.2))+stat_summary(fun.y=median, geom="point", size=5, color="red")+
  scale_y_continuous(limits = c(-0.015, 0.06))
p
ggsave("concordance.pdf")


