rm(list=ls())
## library toolboxes
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(survminer)
library(survival)
##
setwd('C:/Users/ZLiu/Desktop/genefusion/Bioinformatics submission')
##
data1 = read.table('Chimerascan_Patients_groups.csv',header = T, sep=',')
data2 = read.table('SOAPfuse_Patients_groups.csv',header = T, sep=',')
data3 = read.table('Tophat_Patients_groups.csv',header = T, sep=',')
##
data$groups = as.factor(data$groups)
p<-ggplot(data, aes(x=Flag, y=Age, fill=Flag)) +geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2),size=3)+stat_summary(fun.y=median, geom="point", size=2, color="red")
p+scale_color_brewer(palette="Dark2")
