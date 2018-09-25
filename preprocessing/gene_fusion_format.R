rm(list=ls())
# The scirpt is developed to format fusions detected by three 
# calling algorithms including Chimerascan, SOAPfuse and TopHat
# data rearrangement for comparsion
setwd("C:/Users/ZLiu/Desktop/genefusion/Updated_fusion_results_2018.06.05")
data = read.table('Results_Chimerascan_2018.06.05.txt',header=TRUE,sep="\t")
#data = read.table('Results_SOAPfuse_2018.06.05.txt',header=TRUE,sep="\t")
#data = read.table('Results_TopHat_2018.06.05.txt',header=TRUE,sep="\t")
gene1="";gene2="";chr1="";chr2="";patient_ID1="";fusion_index="";
j=0
# handling synonyms in gene1 
for (i in 1:nrow(data)){
  temp = unlist(strsplit(toString(data[i,1]), ","))
  ll=length(temp);
  for (m in 1:ll){
    j=j+1
    gene1[j]=temp[m];gene2[j]=toString(data[i,4]);chr1[j]=toString(data[i,2]);chr2[j]=toString(data[i,5]);patient_ID1[j]=toString(data[i,8]);fusion_index[j]=i
  }
}
rm(temp)
# handling synonyms in gene2
gene11="";gene22="";chr11="";chr22="";patient_ID2="";fusion_index1="";
j=0
for (i in 1:length(gene2)){
  temp = unlist(strsplit(gene2[i], ","))
  ll=length(temp);
  for (m in 1:ll){
    j=j+1
    gene22[j]=temp[m];gene11[j]=gene1[i];chr11[j]=chr1[i];chr22[j]=chr2[i];patient_ID2[j]=patient_ID1[i];fusion_index1[j]=fusion_index[i]
  }
}
rm(gene1,gene2,chr1,chr2,i,j,ll,m,temp,patient_ID1,data)
# sort genes based on alphabet order
gene_fusion_pairs=""
for (i in 1:length(gene11)){
  gene1=paste(gene11[i],'%',chr11[i],sep='')
  gene2=paste(gene22[i],'%',chr22[i],sep='')
  gene_sort=sort(c(gene1,gene2))
  gene_fusion_pairs[i]=paste(gene_sort[1],'%',gene_sort[2],sep='')
}
gene1_Chimerascan=gene11;gene2_Chimerascan=gene22;
gene_fusion_pairs_Chimerascan=gene_fusion_pairs
patient_Chimerascan=patient_ID2
fusion_index_Chimerascan=fusion_index1
save(gene_fusion_pairs_Chimerascan,patient_Chimerascan,fusion_index_Chimerascan,gene1_Chimerascan,gene2_Chimerascan,file = "Chimerascan_fusion.RData")
#gene1_SOAPfuse=gene11;gene2_SOAPfuse=gene22;
#gene_fusion_pairs_SOAPfuse=gene_fusion_pairs
#patient_SOAPfuse=patient_ID2
#fusion_index_SOAPfuse=fusion_index1
#save(gene_fusion_pairs_SOAPfuse,patient_SOAPfuse,fusion_index_SOAPfuse,gene1_SOAPfuse,gene2_SOAPfuse,file = "SOAPfuse_fusion.RData")
#gene1_TopHat=gene11;gene2_TopHat=gene22;
#gene_fusion_pairs_TopHat=gene_fusion_pairs
#patient_TopHat=patient_ID2
#fusion_index_TopHat=fusion_index1
#save(gene_fusion_pairs_TopHat,patient_TopHat,fusion_index_TopHat,gene1_TopHat,gene2_TopHat,file = "TopHat_fusion.RData")
#rm(gene1,gene2,gene_sort,gene11,gene22,chr11,chr22)