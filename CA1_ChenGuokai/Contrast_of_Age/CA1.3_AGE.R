install.packages("pacman")
library(pacman)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
p_load(ggplot2, EBImage, jpeg, ggpubr, plotly)
p_load(pd.mta.1.0, tidyverse, pheatmap, DESeq2, GEOquery, limma, oligo)
library(GEOquery)
library(limma)
library(pd.mta.1.0)
library(oligo)
library(DESeq2)
library(pheatmap)  
library(ggplot2)   
library(tidyverse)

#From Internet
library(AnnoProbe)
p_load(GEOmirror)
remotes::install_github("jmzeng1314/GEOmirror")
library(GEOmirror)


################### PART1::Start to abstract our data  ###################

###In this part, we will get a data frame with gene expression.

geneCELs = list.celfiles('GSE124483_RAW',listGzipped=T,full.name=T)
affyGeneFS <- read.celfiles(geneCELs)   #affyGeneFS:oligoclasses type  #What it means?

genePS <- rma(affyGeneFS, target = "probeset")   ###
featureData(genePS) <- getNetAffx(genePS, "probeset")  ##From Internet Solution

exp1 <- exprs(genePS) #output genePS as matrix
exprdf <- data.frame(exp1)  #transfer matrix into df 

anno=genePS@featureData@data$'geneassignment' #Get gene symbol #why in this geneassignment?

symbol <- str_split_fixed(anno,pattern = "//",3)[,2]  #Cut extra char of gene symbol 
#Usually the second one is the gene name

#To this step, we have already get symbol(name of gene), exprdf(expression of probe)
#Next step is to merge them
temsym<-symbol
temprobe<-exprdf    #temporal data in case modify our parameter

symboldf = data.frame(gene_symbol=c(temsym))  #symboldf describe gene symbol

#I get problem to assign the value of gene name, so i transpose it then to process and it works
ttemprobe=t(temprobe)
colnames(ttemprobe)<-symboldf[,1]
exp10=t(ttemprobe)
#exp10 is the frame that gene alter probe, rownames is gene, colnames is sample
#Next, I did some modification to make it meaningful
exp11=exp10[rownames(exp10)!="",]  #eliminate the probe withour coresponding gene

flag=duplicated(rownames(exp11))#Check repeated gene, use flag as marker, true means repeated
exp12=exp11[flag[]!=T,]   #Select "False" row
####exp111 <- aggregate(.~rownames(exp11),mean,data=exp11)  #试试看能不能根据行名进行
exp13=as.data.frame(exp12) #finally get gene expression data frame exp13
jiyinbiaodaliangbiao<-exp13

#Final, we get our desired data of this part: exp13
#Important variate shouldn't make change later:"exp13"


###################  PART2::Differentially Expression Analysis  ###################

gseineed<-getGEO('GSE124483')#Got troubled with extracting group data from cel file
gseineed<-gseineed[[1]]     #I give up and use getGEO to acquire group design
#pData(gseineed)
#pData(gseineed)[['age:ch1']]
#pData(gseineed)[['genotype/variation:ch1']]
#From annotated code above I read there are 2 most obvious group design taking the age and genotype as factor
#We just use genotype as example to analyse, age type will be the same process

#1.Expression matrix: exprSet
#2.Group matrix: design
#3.Contrast matrix: contrast

# 1.Expression matrix: exprSet
exp<-exp13  

# 2.Group matrix: design matrix
group<-pData(gseineed)[['age:ch1']]   
group <- c(rep("young",3),rep("old",3),rep("young",3),rep("old",3))#Shorten the group name
Group<-as.data.frame(group)
rownames(Group)<-colnames(exp)


designtable<-model.matrix(~0+factor(group))   #[factor] function turn group into facotr 
colnames(designtable)<-levels(factor(group))  #for genotype
rownames(designtable)<-colnames(exp)


# 3.contrast matrix
contrast.matrix=makeContrasts(young - old,levels=designtable)


fit<-lmFit(exp,designtable) 

fit2<-contrasts.fit(fit,contrast.matrix)

fit2<-eBayes(fit2) 


options(digits=4)
DEG<-topTable(fit2,coef="young - old",n=Inf)#coef should be the group name


DEG$group<-ifelse(DEG$P.Value>0.05,"no_change",
                  ifelse(DEG$logFC>1,"up",
                         ifelse(DEG$logFC< -1,"down","no_change")))
#p-value means reliability, so we choose value down to 0.05 group



table(DEG$group)
DEG$gene<-rownames(DEG)

DEG  #Differentially Expressed Gene by genotype

dif <- DEG[DEG[,"P.Value"]<0.05&abs(DEG[,"logFC"])>1,]
guolvhoujuzhen <- dif  #Differentially Expressed Gene (exclude high p-value and unchanged genes)


#Final, we get our desired data of this part: "dif" "DEG" 
#Some important variate shouldn't make change later:"dif" "DEG" "Group" 


###################  PART3:Heatmap  ###################


dif_20<-dif[1:20,]#top20 p.value genes

pzhiqian20 <- dif_20


exp_1<-cbind(exp13,rownames(exp13))  
colnames(exp_1)[13]<-"gene"           #Definte Col13 as gene col
combine<-merge(dif_20,exp_1,by="gene")  #combine is temporal frame. Acodin dif20's gene,find expression value in exp1
rownames(combine)<-combine[,1]       #combine is frame merge generated without rowname #Re-transfer it
com<-combine[,9:20]                 #First 8 col is DEG data, we don't need it anymore
com_n<-as.data.frame(lapply(com, function(x)as.numeric(as.character(x))))#Turn numeric type
rownames(com_n)<-rownames(com)  #last row code also lose the names, rename them


#Start to plot heatmap
df<-com_n
xaxis<-str_split_fixed(colnames(df),pattern = "_",3)[,1] #the file name too long, cut it
colnames(df)<-xaxis  #change cloname
rownames(annotation_col)<-colnames(df)
#annotation_col <- data.frame(group = factor(rep(c("wt", "NLRP3_KO"),each = 6)))


annotation_col <- Group#Group is used again #df's colname correspond Group'srowname, get genetype from Group
#annotation_col <- merge(Group,df,by=intersect()) #abandoned function #搞不懂,不搞了 

###########Cool Cool Heatmapヽ(✿ﾟ▽ﾟ)ノ###########
colnames(df)<-xaxis
retu <- pheatmap(df,
                 display_numbers = T,
                 cluster_rows = T,cluster_cols = F,  #row cluster must be true, otherwises group will be wrong
                 main ="young-old",show_colnames = T,show_rownames = T,
                 legend_breaks = c(7,4,1),legend_labels = c("high","median","low"),
                 fontsize_row = 8,fontsize_col = 10,angle_col = 45,  #angle good, fancy fancy
                 annotation_col = annotation_col
)   
retu


#IGONORE_THIS_PLZ###
#volcanoplot(fit2,coef=1)
#interesting_genes <- topTable(fit2
#                              ,number=Inf,p.value = 0.05,lfc=2)
#
#volcanoplot(fit2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
#points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

###################  PART4:Volcano  ########From Internet###########


need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
huoshantu1 <- deg_volcano(need_deg,1)
huoshantu2 <- deg_volcano(need_deg,2)


#From volcano plot we can find interesting gene.

###################  PART5:Summary  ###################                                   


#important value

#"exprdf" : gene expression data according by probe
#"jiyinbiaodaliangbiao" : gene expression data frame
#"Group" : experiment group design 
#"dif" "guolvhoujuzhen"  : differentially express gene analysis (p<0.05)
#"DEG" "chayibiaodajiyin": differentially express gene analysis (whole)
#"retu" : heatmap of top20 p-value gene
#"huoshantu1" : volcano plot of whole gene(up and down gene noted)
#"huoshantu2" : volcano plot of whole gene(gene name noted)

#If you want to check, you can input above variate to see in console
