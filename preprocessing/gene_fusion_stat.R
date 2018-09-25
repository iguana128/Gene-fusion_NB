rm(list=ls())
# The scirpt is developed to carry out statistical analysis of fusions detected by three 
# calling algorithms including Chimerascan, SOAPfuse and TopHat
# data rearrangement for comparsion

## set directory
setwd("C:/Users/ZLiu/Desktop/genefusion/Updated_fusion_results_2018.06.05")

## load fusion data and annotation data
load('Chimerascan_fusion.RData')
load('SOAPfuse_fusion.RData')
load('TopHat_fusion.RData')
clinical_info = read.table('clinical_information.txt',header=TRUE,sep="\t")
clinical_info = clinical_info[-499,]
protein_Kinases = as.matrix(read.table('human_protein_Kinases_list.txt',header=FALSE,sep="\t",col.names = FALSE))
key_NB_genes = as.matrix(read.table('key_genes_neuroblastoma.txt',header=FALSE,sep="\t",col.names = FALSE))
ChimerDB = read.table('Copy of ChimerDB3.0_ChimerKB.txt',header=TRUE,sep="\t")
TCGA_fusion = read.table('pancanfus.txt',header=TRUE,sep="\t")
# 
gene_fusion_ChimerDB = ""
for (i in 1:nrow(ChimerDB)){
  temp = sort(c(toString(ChimerDB[i,5]),toString(ChimerDB[i,9])))
  gene_fusion_ChimerDB[i] = paste(temp[1],'_',temp[2],sep="")
}
rm(temp)

#
gene_fusion_TCGA = ""
for (i in 1:nrow(TCGA_fusion)){
  temp = sort(c(toString(TCGA_fusion[i,3]),toString(TCGA_fusion[i,4])))
  gene_fusion_TCGA[i] = paste(temp[1],'_',temp[2],sep="")
}
rm(temp)

## Concordance of fusions detected by the three algorithms per patient
patient_id = clinical_info[,1]
tani=matrix(0,nrow(clinical_info),1)
for (i in 1:nrow(clinical_info)){
  v1 = gene_fusion_pairs_Chimerascan[which(patient_Chimerascan==patient_id[i])]
  v2 = gene_fusion_pairs_SOAPfuse[which(patient_SOAPfuse==patient_id[i])]
  v3 = gene_fusion_pairs_TopHat[which(patient_TopHat==patient_id[i])]
  A=setdiff(v1,union(v2,v3));B=setdiff(v2,union(v1,v3));C=setdiff(v3,union(v1,v2))
  AB=setdiff(intersect(v1,v2),v3);BC=setdiff(intersect(v2,v3),v1);AC=setdiff(intersect(v1,v3),v2)
  ABC=intersect(v1,intersect(v2,v3))
  tani[i]=(length(AB)+length(AC)+length(BC)+length(ABC))/(length(A)+length(B)+length(C)+length(AB)+length(AC)+length(BC)+length(ABC))
  rm(v1,v2,v3,A,B,C,AB,AC,BC,ABC)
}

## assign the clinical information into patients in different fusion algorithms
clinical_Chimerascan_index = matrix("",length(patient_Chimerascan),4)
clinical_SOAPfuse_index = matrix("",length(patient_SOAPfuse),4)
clinical_TopHat_index = matrix("",length(patient_TopHat),4)
#
hh=match(patient_Chimerascan,clinical_info[,1])
clinical_Chimerascan_index[,1]=as.character(clinical_info[hh,3])
clinical_Chimerascan_index[,2]=as.character(clinical_info[hh,4])
clinical_Chimerascan_index[,3]=as.character(clinical_info[hh,5])
clinical_Chimerascan_index[,4]=as.character(clinical_info[hh,6])
rm(hh)
#
hh=match(patient_SOAPfuse,clinical_info[,1])
clinical_SOAPfuse_index[,1]=as.character(clinical_info[hh,3])
clinical_SOAPfuse_index[,2]=as.character(clinical_info[hh,4])
clinical_SOAPfuse_index[,3]=as.character(clinical_info[hh,5])
clinical_SOAPfuse_index[,4]=as.character(clinical_info[hh,6])
rm(hh)
#
hh=match(patient_TopHat,clinical_info[,1])
clinical_TopHat_index[,1]=as.character(clinical_info[hh,3])
clinical_TopHat_index[,2]=as.character(clinical_info[hh,4])
clinical_TopHat_index[,3]=as.character(clinical_info[hh,5])
clinical_TopHat_index[,4]=as.character(clinical_info[hh,6])
rm(hh)

## annotate the fusion based on protein Kinase family cureated from uiport website
Kinase_fusion_Chimerascan_index = matrix(0,length(patient_Chimerascan),1)
Kinase_fusion_SOAPfuse_index = matrix(0,length(patient_SOAPfuse),1)
Kinase_fusion_TopHat_index = matrix(0,length(patient_TopHat),1)

Kinase_fusion_Chimerascan = union(which(gene1_Chimerascan %in% intersect(gene1_Chimerascan,protein_Kinases)),which(gene2_Chimerascan %in% intersect(gene2_Chimerascan,protein_Kinases)))
Kinase_fusion_SOAPfuse = union(which(gene1_SOAPfuse %in% intersect(gene1_SOAPfuse,protein_Kinases)),which(gene2_SOAPfuse %in% intersect(gene2_SOAPfuse,protein_Kinases)))
Kinase_fusion_TopHat = union(which(gene1_TopHat %in% intersect(gene1_TopHat,protein_Kinases)),which(gene2_TopHat %in% intersect(gene2_TopHat,protein_Kinases)))

Kinase_fusion_Chimerascan_index[Kinase_fusion_Chimerascan] = 1
Kinase_fusion_SOAPfuse_index[Kinase_fusion_SOAPfuse] = 1
Kinase_fusion_TopHat_index[Kinase_fusion_TopHat] = 1

rm(Kinase_fusion_Chimerascan,Kinase_fusion_SOAPfuse,Kinase_fusion_TopHat)

## annotate the fusion based on reported by fusion-related knowledge base (ChimerDB 3.0)
report_fusion_Chimerascan_index = matrix(0,length(patient_Chimerascan),1)
report_fusion_SOAPfuse_index = matrix(0,length(patient_SOAPfuse),1)
report_fusion_TopHat_index = matrix(0,length(patient_TopHat),1)

reported_fusion_Chimerscan_ChimerDB = which(gene_fusion_Chimerscan %in% intersect(gene_fusion_ChimerDB,gene_fusion_Chimerscan))
reported_fusion_SOAPfuse_ChimerDB = which(gene_fusion_SOAPfuse %in% intersect(gene_fusion_ChimerDB,gene_fusion_SOAPfuse))
reported_fusion_TopHat_ChimerDB = which(gene_fusion_TopHat %in% intersect(gene_fusion_ChimerDB,gene_fusion_TopHat))

reported_fusion_Chimerscan_TCGA = which(gene_fusion_Chimerscan %in% intersect(gene_fusion_TCGA,gene_fusion_Chimerscan))
reported_fusion_SOAPfuse_TCGA = which(gene_fusion_SOAPfuse %in% intersect(gene_fusion_TCGA,gene_fusion_SOAPfuse))
reported_fusion_TopHat_TCGA = which(gene_fusion_TopHat %in% intersect(gene_fusion_TCGA,gene_fusion_TopHat))

report_fusion_Chimerascan_index[union(reported_fusion_Chimerscan_ChimerDB,reported_fusion_Chimerscan_TCGA)] = 1
report_fusion_SOAPfuse_index[union(reported_fusion_SOAPfuse_ChimerDB,reported_fusion_SOAPfuse_TCGA)] = 1
report_fusion_TopHat_index[union(reported_fusion_TopHat_ChimerDB,reported_fusion_TopHat_TCGA)] = 1
rm(reported_fusion_Chimerscan_TCGA,reported_fusion_SOAPfuse_TCGA,reported_fusion_TopHat_TCGA)
rm(reported_fusion_Chimerscan_ChimerDB,reported_fusion_SOAPfuse_ChimerDB,reported_fusion_TopHat_ChimerDB)

## annotate the fusion by key genes in neuroblastoma curated by literature survery
key_NB_genes_fusion_Chimerascan_index = matrix(0,length(patient_Chimerascan),1)
key_NB_genes_fusion_SOAPfuse_index = matrix(0,length(patient_SOAPfuse),1)
key_NB_genes_fusion_TopHat_index = matrix(0,length(patient_TopHat),1)

key_NB_genes_fusion_Chimerascan = union(which(gene1_Chimerascan %in% intersect(gene1_Chimerascan,key_NB_genes)),which(gene2_Chimerascan %in% intersect(gene2_Chimerascan,key_NB_genes)))
key_NB_genes_fusion_SOAPfuse = union(which(gene1_SOAPfuse %in% intersect(gene1_SOAPfuse,key_NB_genes)),which(gene2_SOAPfuse %in% intersect(gene2_SOAPfuse,key_NB_genes)))
key_NB_genes_fusion_TopHat = union(which(gene1_TopHat %in% intersect(gene1_TopHat,key_NB_genes)),which(gene2_TopHat %in% intersect(gene2_TopHat,key_NB_genes)))

key_NB_genes_fusion_Chimerascan_index[key_NB_genes_fusion_Chimerascan] = 1
key_NB_genes_fusion_SOAPfuse_index[key_NB_genes_fusion_SOAPfuse] = 1
key_NB_genes_fusion_TopHat_index[key_NB_genes_fusion_TopHat] = 1
rm(key_NB_genes_fusion_Chimerascan,key_NB_genes_fusion_SOAPfuse,key_NB_genes_fusion_TopHat)

## output result table
#
temp =matrix("",length(gene_fusion_pairs_Chimerascan),4)
for (i in 1:length(gene_fusion_pairs_Chimerascan)){
  temp[i,] = unlist(strsplit(gene_fusion_pairs_Chimerascan[i], "%"))
}
Chimerascan_results=cbind(patient_Chimerascan,gene_fusion_Chimerscan,temp,Kinase_fusion_Chimerascan_index,report_fusion_Chimerascan_index,key_NB_genes_fusion_Chimerascan_index,clinical_Chimerascan_index)
colnames(Chimerascan_results)=c('patient_Chimerascan','gene_fusion_Chimerscan','gene1','chr1','gene2','chr2','Kinase_related_fusion','reported_fusion','key_NB_fusion','INSS_stage','AGE','High_risk','MYCN_status')
save(Chimerascan_results,file='Chimerascan_results.RData')
rm(temp)
#
temp =matrix("",length(gene_fusion_pairs_SOAPfuse),4)
for (i in 1:length(gene_fusion_pairs_SOAPfuse)){
  temp[i,] = unlist(strsplit(gene_fusion_pairs_SOAPfuse[i], "%"))
}
SOAPfuse_results=cbind(patient_SOAPfuse,gene_fusion_SOAPfuse,temp,Kinase_fusion_SOAPfuse_index,report_fusion_SOAPfuse_index,key_NB_genes_fusion_SOAPfuse_index,clinical_SOAPfuse_index)
colnames(SOAPfuse_results)=c('patient_SOAPfuse','gene_fusion_SOAPfuse', 'gene1','chr1','gene2','chr2','Kinase_related_fusion','reported_fusion','key_NB_fusion','INSS_stage','AGE','High_risk','MYCN_status')
save(SOAPfuse_results,file='SOAPfuse_results.RData')
rm(temp)
#
temp =matrix("",length(gene_fusion_pairs_TopHat),4)
for (i in 1:length(gene_fusion_pairs_TopHat)){
  temp[i,] = unlist(strsplit(gene_fusion_pairs_TopHat[i], "%"))
}
TopHat_results=cbind(patient_TopHat,gene_fusion_TopHat,temp,Kinase_fusion_TopHat_index,report_fusion_TopHat_index,key_NB_genes_fusion_TopHat_index,clinical_TopHat_index)
colnames(TopHat_results)=c('patient_TopHat','gene_fusion_TopHat','gene1','chr1','gene2','chr2','Kinase_related_fusion','reported_fusion','key_NB_fusion','INSS_stage','AGE','High_risk','MYCN_status')
save(TopHat_results,clinical_info,tani,file='TopHat_results.RData')
rm(temp)






