rm(list=ls())
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(survminer)
library(survival)
##
setwd("C:/Users/ZLiu/Desktop/genefusion/Final_results_forHCA/Final_results_forHCA")
##
data = read.table('Chimerascan_Patients_groups.csv',header = T, sep=',')
data$groups = as.factor(data$groups)
p<-ggplot(data, aes(x=groups, y=OS_d, fill=groups)) +geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2),size=3)+stat_summary(fun.y=median, geom="point", size=2, color="red")
p+scale_color_brewer(palette="Dark2")
ggsave('figure3a1.pdf')
rm(data)
##
data = read.table('SOAPfuse_Patients_groups.csv',header = T, sep=',')
data$groups = as.factor(data$groups)
p<-ggplot(data, aes(x=groups, y=OS_d, fill=groups)) +geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2),size=3)+stat_summary(fun.y=median, geom="point", size=2, color="red")
p+scale_color_brewer(palette="Dark2")
ggsave('figure3a2.pdf')
rm(data)
##
data = read.table('Tophat_Patients_groups.csv',header = T, sep=',')
data$groups = as.factor(data$groups)
p<-ggplot(data, aes(x=groups, y=OS_d, fill=groups)) +geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2),size=3)+stat_summary(fun.y=median, geom="point", size=2, color="red")
p+scale_color_brewer(palette="Dark2")
ggsave('figure3a3.pdf')
rm(data)
##
data = read.table('Chimerascan_patient_relation_matrix.csv',header = T, sep=',')
data1=data[,-1]
row.names(data1) = paste(rep('V',length(data[,1])),data[,1],sep='')
melt_data = melt(as.matrix(data1))
mypalette<- colorRampPalette(c("white", "darkblue"))(n = 10)
heatmap(as.matrix(data1),hclustfun = hclust,col = mypalette)
rm(data,data1,melt_data)
##
data = read.table('SOAPfuse_patient_relation_matrix.csv',header = T, sep=',')
data1=data[,-1]
row.names(data1) = paste(rep('V',length(data[,1])),data[,1],sep='')
melt_data = melt(as.matrix(data1))
mypalette<- colorRampPalette(c("white", "darkblue"))(n = 5)
heatmap(as.matrix(data1),hclustfun = hclust,col = mypalette)
rm(data,data1,melt_data)
##
data = read.table('TopHat_patient_relation_matrix.csv',header = T, sep=',')
data1=data[,-1]
row.names(data1) = paste(rep('V',length(data[,1])),data[,1],sep='')
melt_data = melt(as.matrix(data1))
mypalette<- colorRampPalette(c("white", "darkblue"))(n = 5)
heatmap(as.matrix(data1),hclustfun = hclust,col = mypalette)
rm(data,data1,melt_data)

## survival analysis
setwd("C:/Users/ZLiu/Desktop/genefusion/Updated_results_2018.06.13")
data = read.table('ChimeraScan_survival.csv',header = T, sep=',')
fit <- survfit(Surv(OS_days, OS_bin) ~ Flag, data = data)
ggsurvplot(fit, data = data, risk.table = TRUE)
#
data = read.table('SOAPfuse_survival.csv',header = T, sep=',')
fit <- survfit(Surv(OS_days, OS_bin) ~ Flag, data = data)
ggsurvplot(fit, data = data, risk.table = TRUE)
#
data = read.table('TopHat_survival.csv',header = T, sep=',')
fit <- survfit(Surv(OS_days, OS_bin) ~ Flag, data = data)
ggsurvplot(fit, data = data, risk.table = TRUE)



