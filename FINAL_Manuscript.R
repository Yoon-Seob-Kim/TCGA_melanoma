##Library Load
library(devtools)
library(factoextra)
library(SigProfilerExtractor)
library(edgeR)
library(ggpubr)
library(maftools)
library(data.table)
library(ggbiplot)
library(NMF)
library(survival)
library(survminer)
library(dplyr)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(easyGgplot2)
library(ggrepel)
library(pheatmap)
library(tidyr)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(maftools)
library(dndscv)
library(GSVA)
library(edgeR)
library(ggrepel)
library(ggpubr)
library(copynumber)
library(rngtools)
library(caTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(deconstructSigs)
library(MutationalPatterns)

##visualization 96-mutational context of each sample 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/00_RAW/TCGA")
sample=read.delim("list.txt",sep="\t",header=F)
vcf_files =sample$V1
sample_names=sample$V2
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
fwrite(x=as.data.frame(mut_mat),file="context.csv",row.names=T,col.names=T)

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/00_RAW/ICGC_WES")
sample=read.delim("list.txt",sep="\t",header=F)
vcf_files =sample$V1
sample_names=sample$V2
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
fwrite(x=as.data.frame(mut_mat),file="context.csv",row.names=T,col.names=T)


setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/00_RAW/SRA")
sample=read.delim("list.txt",sep="\t",header=F)
vcf_files =sample$V1
sample_names=sample$V2
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
fwrite(x=as.data.frame(mut_mat),file="context.csv",row.names=T,col.names=T)


##Run sigprofiler
#Figure S1
sigprofilerextractor("matrix", "/data/Delta_data4/kysbbubbu/Frontier_TCGA/01_Sigprofiler/TCGA", "/data/Delta_data4/kysbbubbu/Frontier_TCGA/01_Sigprofiler/Sigprofiler.txt",stability=0.9, exome=T,context_type="96",cosmic_version=3.2,maximum_signatures=8)

##Perform decomposition 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_NMF")
NMF=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
NMF=t(NMF)
res=nmf(NMF, 2, nrun=1000, seed=123456, .options="p24")
#Export basis and coefficient
basis=as.matrix(res@fit@W)
basis=t(t(basis)/colSums(basis))
write.csv(x=basis, file="basis.csv", row.names=TRUE)
coefmap=as.matrix(res@fit@H)
coefmap=t(t(coefmap)/colSums(coefmap))
write.csv(x=coefmap, file="coefmap.csv", row.names=TRUE)
#K-means clustering 
cluster=kmeans(t(coefmap)[,1:2],2)$cluster
cluster=as.factor(cluster)
current.cluster.ids <- c("1","2")
new.cluster.ids <- c("UV-high","UV-low")
cluster  <- plyr::mapvalues(x =cluster, from = current.cluster.ids, to = new.cluster.ids)
cluster=factor(cluster, levels= c("UV-low","UV-high"))
summary(cluster)
write.csv(x=as.data.frame(cluster), file="cluster.csv", row.names=TRUE)

#Define cosine similarity
cos.sim=function(ma, mb){
  mat=tcrossprod(ma, mb)
  t1=sqrt(apply(ma, 1, crossprod))
  t2=sqrt(apply(mb, 1, crossprod))
  mat / outer(t1,t2)}

#perform deconstrucsig
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig = read.table("deconstruct.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig=as.data.frame(t(sig))

#Figure 1. De novo mutational signature extraction and unsupervised clustering in TCGA CM dataset. 
##Figure 1A 
SigA=whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, sample.id = 'SigA',associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
cos.sim(SigA$tumor,SigA$product)##cos 0.952
pdf('SigA_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SigA)
dev.off()

SigB=whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019, contexts.needed = FALSE, sample.id = 'SigB',associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
cos.sim(SigB$tumor,SigB$product)##cos 0.999
pdf('SigB_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SigB)
dev.off()

#comparison between COSMIC SBS signatures and NMF-decomposed basis
sig=as.data.frame(t(basis))
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)
COSMIC.v3=signatures.exome.cosmic.v3.may2019

#Cosine similarity between SigA and COSMIC signatures
SigA_cos=as.matrix(sig[1,])
COSMIC.v3=as.matrix(signatures.exome.cosmic.v3.may2019)
cos=cos.sim(SigA_cos,COSMIC.v3)
cos=round(cos[,order(cos,decreasing=TRUE)])
cos

#Cosine similarity between SigB and COSMIC signatures
SigB_cos=as.matrix(sig[2,])
COSMIC.v3=as.matrix(signatures.exome.cosmic.v3.may2019)
cos2=cos.sim(SigB_cos,COSMIC.v3)
cos2=round(cos2[,order(cos2,decreasing=TRUE)],2)
cos2

##Piechart 
barplot = as.data.frame(rep(c("SBS1", "SBS3", "SBS5", "SBS7a", "SBS7b","SBS38","Unkown"),2))
colnames(barplot)="SBS"
barplot$proportion=c("0.205","0.321","0.222","0","0","0.089","0.163","0","0","0.076","0.512","0.405","0","0.007")
barplot$basis=rep(c("SigA", "SigB"),each=7)
barplot$proportion=as.numeric(as.character(barplot$proportion))
barplot$SBS=factor(barplot$SBS, levels= c("SBS1","SBS3","SBS5","SBS7a","SBS7b","SBS38","Unkown"))
#Figure 1B
ggplot(data = barplot, aes(x = "", y = proportion, fill = SBS )) +   geom_bar(stat = "identity", position = position_fill()) +   theme_minimal()+ theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold")) + coord_polar(theta = "y") +  facet_wrap(~ basis)  + theme(axis.text=element_blank(),legend.position='top',legend.title=element_blank(),legend.text=element_text(size=12))

###Figure 1C
#Heatmap showing proportion of SigA/SigB
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
coefmap2 = coefmap[,order(dat.n$SigB)]
rownames(coefmap2)=c("SigA","SigB")
my_gene_col <- data.frame(sample = rep(c("SigA", "SigB"), c(1,1)))
my_sample_col=as.data.frame(dat.n$cluster)
colnames(my_sample_col)="Cluster"
rownames(my_sample_col)=colnames(coefmap)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels= c("UV-low", "UV-high"))
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Cluster=Var2)
pheatmap(coefmap2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_names_row=FALSE,annotation_col =my_sample_col,annotation_names_col=FALSE,annotation_colors=ann_colors)

#Bar plot of log-scale mutation count 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
colnames(dat.n)
dat.n2=dat.n[order(dat.n$SigB),]
ggplot(dat.n2, aes(x=factor(1:length(Mutation)), y=Mutation)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

#Heatmap showing 96-class mutational context
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n2=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="Cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$Cluster=dat.n2$cluster
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
dat.n = dat.n[,order(dat.n2$SigB)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Type = Var1, Cluster=Var2)
##900/600
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE)


##Supplementary Figure S2. Average mutational context of total TCGA cohort, UV-low and UV-high clusters. 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig = read.table("deconstruct.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig=as.data.frame(t(sig))
Total = whichSignatures(tumor.ref = sig, signatures.ref= signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, sample.id = 'Total', associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
cos.sim(Total$tumor,Total$product)##cos 0.996
pdf('Total_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(Total)
dev.off()

UV_low = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'UV_low')
cos.sim(UV_low$tumor,UV_low$product)##cos 0.954
pdf('UV_low_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(UV_low)
dev.off()

UV_high = whichSignatures(tumor.ref = sig, signatures.ref= signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, sample.id = 'UV_high', associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
cos.sim(UV_high$tumor,UV_high$product)##cos 0.999
pdf('UV_high_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(UV_high)
dev.off()

##Supplementary Figure S3. Comparison of mutational context using variants called by Mutect2 and Varscan2. 
#Mutational context from Mutect2 data 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n2=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n= dat.n[,order(dat.n2$Heatmap_order)]
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="Cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$Cluster=rep(c("UV-low","UV-high"), c(75,391))
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend	=FALSE)

#Mutational context from Varscan2 data
dat.n=read.table("Varscan2_heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n2=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n=dat.n[,order(dat.n2$Heatmap_order)]
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$cluster=rep(c("UV-low","UV-high"), c(75,391))
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="type"
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(type = Var1, cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend=TRUE)

#Cosine similarity between Mutect2 and Varscan2  
dat.n=as.matrix(read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE))
dat.n2=as.matrix(read.table("Varscan2_heatmap.txt", sep="\t", row.names=1, header=TRUE))
cos=cos.sim(dat.n,dat.n2)
cos.matrix=as.data.frame(diag(cos))
dat.n2=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
cos.matrix2=as.data.frame(cos.matrix[order(dat.n2$Heatmap_order),])
colnames(cos.matrix2)="cos.sim"
ggplot(cos.matrix2, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

##Supplementary Figure S4. Comparison of average read depth, proportion of neoadjuvant chemotherapy, and primary tumor diagnosis between UV-low and UV-high clusters. 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
my_comparisons <- list( c("UV-low","UV-high"))
ggboxplot(dat.n, x = "cluster", y = "Read_depth.1",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text= element_text(size = 15,face="bold",color="black"), axis.title.y = element_text(size=17, face="bold"))+theme(legend.title =element_text(size = 16,face="bold"),legend.text = element_text(size=16), legend.position="right")+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +ylab("Read depth (ENSEMBL exonic region)") +labs(colour="Cluster")+ theme(legend.position = "none")
wilcox.test(Read_depth.1 ~ cluster, data = dat.n)

dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
dat.n2 <- as.data.frame(table(dat.n$cluster,dat.n$Neoadjuvant.CTx.Hx))
dat.n2$Var2=factor(dat.n2$Var2,levels=c("No","Yes"))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text= element_text(size = 15,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title =element_text(size = 16,face="bold"),legend.text = element_text(size=16), legend.position="right")+labs(fill="Neoadjuvant CTx")
dat.n2
chisq.test(matrix(c(6,41,69,350),nrow = 2, ncol = 2))

dat.n2 <- as.data.frame(table(dat.n$cluster,dat.n$Primary.tumor.diagnosed))
dat.n2$Var2=factor(dat.n2$Var2,levels=c("No","Yes"))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text= element_text(size = 15,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title =element_text(size = 16,face="bold"),legend.text = element_text(size=16), legend.position="right")+labs(fill="Primary tumor")
chisq.test(matrix(c(69,372,6,19),nrow = 2, ncol = 2))

###NMF decomposition of SRA dataset
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
SRA_NMF=read.table("SRA_heatmap.txt", sep="\t", row.names=1, header=TRUE)
SRA_NMF = t(SRA_NMF)
SRA_res <- nmf(SRA_NMF, 2, nrun=1000, seed=123456, .options="p24")
#Export basis and coefficient
basis=as.matrix(SRA_res@fit@W)
basis=t(t(basis)/colSums(basis))
write.csv(x=basis, file="SRA_basis.csv", row.names=TRUE)
SRA_coefmap=as.matrix(SRA_res@fit@H)
SRA_coefmap=t(t(SRA_coefmap)/colSums(SRA_coefmap))
write.csv(x=SRA_coefmap, file="SRA_coefmap.csv",row.names=TRUE)

###NMF decomposition of ICGC dataset
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
NMF2=read.table("ICGC_heatmap.txt", sep="\t", row.names=1, header=TRUE)
NMF2 = t(NMF2)
res2 <- nmf(NMF2, 2, nrun=1000, seed=123456, .options="p24")
#Export basis and coefficient
basis2=as.matrix(res2@fit@W)
basis2=t(t(basis2)/colSums(basis2))
write.csv(x=basis2, file="ICGC_basis.csv",row.names=TRUE)
icgc_coefmap=as.matrix(res2@fit@H)
icgc_coefmap=t(t(icgc_coefmap)/colSums(icgc_coefmap))
write.csv(x=icgc_coefmap, file="ICGC_coefmap.csv",row.names=TRUE)

#Comparison between ICGC basis and TCGA basis
basis2=as.matrix(res2@fit@W)
basis2=t(t(basis2)/colSums(basis2))
basis2=t(basis2)
colnames(basis2)=colnames(signatures.exome.cosmic.v3.may2019)
sig=as.matrix(res@fit@W)
sig=t(t(sig)/colSums(sig))
sig=t(sig)
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)
ICGC_Basis1_cos=t(as.matrix(basis2[1,]))
sig_cos=t(as.matrix(sig[1,]))
cos.sim(ICGC_Basis1_cos,sig_cos)
ICGC_Basis2_cos=t(as.matrix(basis2[2,]))
sig_cos=t(as.matrix(sig[2,]))
cos.sim(ICGC_Basis2_cos,sig_cos)


#SRA basis --TCGA basis comparison 
basis=as.matrix(res@fit@W)
basis=t(t(basis)/colSums(basis))
sig=t(basis)
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)

basis=as.matrix(SRA_res@fit@W)
basis=t(t(basis)/colSums(basis))
SRA_sig=t(basis)
colnames(SRA_sig)=colnames(signatures.exome.cosmic.v3.may2019)

SRA_Basis1_cos=t(as.matrix(SRA_sig[1,]))
Basis1_cos=t(as.matrix(sig[1,]))
cos.sim(SRA_Basis1_cos,Basis1_cos)

SRA_Basis2_cos=t(as.matrix(SRA_sig[2,]))
Basis2_cos=t(as.matrix(sig[2,]))
cos.sim(SRA_Basis2_cos,Basis2_cos)

##Supplementary Figure S5. NMF-based unsupervised clustering of mutational signatures in ICGC and SRA datasets. 
#Figure S5a. ICGC basis visulization
basis=as.matrix(res2@fit@W)
basis=t(t(basis)/colSums(basis))
icgc_Basis=as.data.frame(t(basis))
colnames(icgc_Basis)=colnames(signatures.exome.cosmic.v3.may2019)
rownames(icgc_Basis)=c("Basis1","Basis2")

##cos 0.955
icgc_Basis1 = whichSignatures(tumor.ref = icgc_Basis, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'Basis1')
cos.sim(icgc_Basis1$tumor,icgc_Basis1$product)

##cos 0.998
icgc_Basis2 = whichSignatures(tumor.ref = icgc_Basis, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'Basis2')
cos.sim(icgc_Basis2$tumor,icgc_Basis2$product)


pdf('icgc_Basis1_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(icgc_Basis1)
dev.off()

pdf('icgc_Basis2_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(icgc_Basis2)
dev.off()

#Figure S5b. SRA basis visulization
SRA_Basis=as.data.frame(rbind(SRA_Basis1_cos,SRA_Basis2_cos))
rownames(SRA_Basis)=c("Basis1","Basis2")
colnames(SRA_Basis)=colnames(sigs.input)

##cos 0.998
SRA_Basis1 = whichSignatures(tumor.ref = SRA_Basis, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'Basis1')
cos.sim(SRA_Basis1$tumor,SRA_Basis1$product)

##cos 0.942
SRA_Basis2 = whichSignatures(tumor.ref = SRA_Basis, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'Basis2')
cos.sim(SRA_Basis2$tumor,SRA_Basis2$product)

pdf('SRA_Basis1_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SRA_Basis1)
dev.off()

pdf('SRA_Basis2_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SRA_Basis2)
dev.off()

##Figure S5c. ICGC dataset coef/mutation/context visualization
icgc_annot=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE, quote = "")
icgc_annot=icgc_annot[order(icgc_annot$Basis2_NMF),]
my_sample_col=as.data.frame(icgc_annot[,11])
icgc_annot=icgc_annot[,9:10]
dat.n=t(icgc_annot)
rownames(dat.n)=c("SigA","SigB")
colnames(my_sample_col)="Cluster"
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Cluster=Var2)
rownames(my_sample_col)=colnames(dat.n)
##1000/150
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,annotation_legend = FALSE,legend=FALSE)

icgc_annt=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
icgc_annt=icgc_annt[order(icgc_annt$Basis2_NMF),]
ggplot(icgc_annt, aes(x=factor(1:length(WES_SNP)), y=WES_SNP)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

dat.n=read.table("ICGC_heatmap.txt", sep="\t", row.names=1, header=TRUE, quote = "")
dat.n=t(dat.n)
my_gene_col <- data.frame(Type = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
SRA_annt=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
colnames(SRA_annt)
my_sample_col=as.data.frame(SRA_annt[,"dataset"])
head(SRA_annt)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Study"
dat.n2=dat.n[,order(SRA_annt$Basis2_NMF)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 =c("#E41A1C", "#377EB8")
names(Var2) = c("SKCM-US", "MELA-AU")
ann_colors = list(Type = Var1, Study=Var2)
my_sample_col$Study
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col=my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)

dat.n=read.table("ICGC_heatmap.txt", sep="\t", row.names=1, header=TRUE, quote = "")
icgc_annt=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(icgc_annt$Basis2_NMF),]
sum=as.data.frame(cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96])))
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% gather(variable, value, -ind)
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")

##Figure S5d. SRA dataset coef/mutation/context visualization
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
SRA_annt=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
SRA_annt2=SRA_annt[order(SRA_annt$Basis2_NMF),]
SRA_annt3=t(SRA_annt2[,8:9])
my_sample_col=as.data.frame(SRA_annt[,10])
rownames(my_sample_col)=rownames(SRA_annt)
colnames(my_sample_col)="Cluster"
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Cluster=Var2)
pheatmap(SRA_annt3,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend=FALSE,annotation_colors=ann_colors,legend=FALSE)

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/05_VALIDATION")
SRA_annt=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
SRA_annt2=SRA_annt[order(SRA_annt$Basis2_NMF),]
ggplot(SRA_annt2, aes(x=factor(1:length(Mutation)), y=Mutation)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

dat.n=as.data.frame(t(read.table("SRA_heatmap.txt", sep="\t", row.names=1, header=TRUE)))
my_gene_col <- data.frame(Type = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
SRA_annt=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col=as.data.frame(SRA_annt[,3])
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Study"
dat.n=dat.n[,order(SRA_annt$Basis2_NMF)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = brewer.pal(7,"Set1")
names(Var2) = c("Hugo", "Krauthammer", "Liang", "Riaz", "Roh", "Shain", "Snyder")
ann_colors = list(Type = Var1, Study=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col=my_sample_col,color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors)

dat.n=read.table("SRA_heatmap.txt", sep="\t", row.names=1, header=TRUE)
SRA_annt=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(SRA_annt$Basis2_NMF),]
sum=as.data.frame(cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96])))
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% gather(variable, value, -ind)
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")


##Supplementary Figure S6. Comparison on mutation counts in each dataset (TCGA, ICGC, and SRA). 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
SRA_annt=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
SRA_annt$cluster=factor(SRA_annt$cluster, levels=c("UV-low","UV-high"))
my_comparisons <- list( c("UV-low","UV-high") )
ggboxplot(SRA_annt, x = "cluster", y = "Mutation",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 12,face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
wilcox.test(Mutation ~ NMF_cluster, data = SRA_annt)

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_VALIDATION")
SRA_annt=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
colnames(SRA_annt)
SRA_annt$cluster=factor(SRA_annt$cluster, levels=c("UV-low","UV-high"))
ggboxplot(SRA_annt, x = "cluster", y = "Mutation",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text( size =  16,face="bold",color="black"),axis.text.y = element_text( size = 12,face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
wilcox.test(Mutation ~ cluster, data = SRA_annt)

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_VALIDATION")
SRA_annt=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
colnames(SRA_annt)
SRA_annt$cluster=factor(SRA_annt$cluster, levels=c("UV-low","UV-high"))
ggboxplot(SRA_annt, x = "cluster", y = "Mutation",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text( size =  16,face="bold",color="black"),axis.text.y = element_text( size = 12,face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
wilcox.test(Mutation ~ cluster, data = SRA_annt)

##Supplementary Figure S7. Principal component analysis plot of mutational context of three datasets (TCGA, ICGC, and SRA). 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n=as.matrix(dat.n)
icgc_annt=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n_pca <- prcomp(dat.n, center = TRUE,scale. = FALSE)
fviz_pca_ind(dat.n_pca, col.ind = factor(icgc_annt$cluster),palette = c("red",  "blue"),label="none",title="")+ theme_bw()+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ theme(legend.position = "none")

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
dat.n=read.table("ICGC_heatmap.txt", sep="\t", row.names=1, header=TRUE, quote = "")
icgc_annt=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(icgc_annt$Basis2_NMF),]
dat.n = as.matrix(dat.n)
dat.n_pca <- prcomp(dat.n, center = TRUE,scale. = FALSE)
fviz_pca_ind(dat.n_pca, col.ind = factor(icgc_annt$NMF_cluster),palette = c("red",  "blue"),label="none",title="")+ theme_bw()+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ theme(legend.position = "none")

SRA_NMF=read.table("SRA_heatmap.txt", sep="\t", row.names=1, header=TRUE)
SRA_NMF = as.matrix(SRA_NMF)
SRA_pca <- prcomp(SRA_NMF, center = TRUE,scale. = FALSE)
SRA_annt=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
fviz_pca_ind(SRA_pca,       col.ind = factor(SRA_annt$NMF_cluster),palette = c("red",  "blue"),label="none",title="")+ theme_bw()+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ theme(legend.position = "none")

##Supplementary Figure S8. Comparison of mutational context of all variants and variants within exonic region. 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
dat.n=read.table("ICGC_heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=read.table("ICGC_WGS_heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n=dat.n[1:198,]
dat.n2=dat.n2[rownames(dat.n), ]
dat.n3=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n3=dat.n3[1:198,]
dat.n3=dat.n3[order(dat.n3$Basis2_NMF),]
dat.n=dat.n[rownames(dat.n3), ]
dat.n2=dat.n2[rownames(dat.n3), ]
dat.n=as.matrix(dat.n)
dat.n2=as.matrix(dat.n2)
cos=cos.sim(dat.n,dat.n2)
cos.matrix=diag(cos)
cos.matrix=as.data.frame(cos.matrix)
min(cos.matrix)
max(cos.matrix)
colnames(cos.matrix)="cos.sim"
median(cos.matrix$cos.sim)
ggplot(cos.matrix, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

proportion=as.data.frame(dat.n3$WES_SNP/dat.n3$WGS_SNP*100)
colnames(proportion)="proportion"
ggplot(proportion, aes(x=factor(1:length(proportion)), y=proportion)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Proportion of SNVs")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))
median(proportion$proportion)

my_sample_col=as.data.frame(dat.n3$NMF_cluster)
rownames(my_sample_col)=rownames(dat.n3)
colnames(my_sample_col)="Cluster"
dat.n=t(dat.n)
my_gene_col <- data.frame(Type = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_row =my_gene_col,annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend=FALSE,annotation_colors=ann_colors,legend=FALSE)
dat.n2=t(dat.n2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend=FALSE,annotation_colors=ann_colors,legend=FALSE)

###TCGA cohort deconstrucsig (COSMIC signature refitting)
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
sig=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
sig=as.data.frame(sig)
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)
list=rownames(sig)
for (list1 in list) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct/decon_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct/decon_sig")
  weight=dec$weights
  write.table(x=weight, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

###TCGA cohort deconstrucsig (Basis signature refitting)
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig2 = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(signatures.exome.cosmic.v3.may2019)
sig=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
sig=as.data.frame(sig)
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)
list=rownames(sig)
for (list1 in list) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =sig2, contexts.needed = FALSE, associated	=c("Basis1", "Basis2"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct/basis_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct/basis_sig")
  weight=dec$weights
  write.table(x=weight, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}


###Figure 2. Signature refitting analysis. 
#Figure 2A
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
dat.n=read.table("COSMIC_signature.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n=dat.n[1:20,]
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$cluster=rep(c("UV-low", "UV-high"), c(75,391))
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE,fontsize = 12)
?pheatmap

dat.n=read.table("COSMIC_signature.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n=dat.n[26:28,]
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$cluster=rep(c("UV-low", "UV-high"), c(75,391))
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
rownames(dat.n)=c("SigA","SigB","Unknown")
ann_colors = list(cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors,legend=F,annotation_legend=F)

####TCGA Age distribution 
####Figure2B
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$age=as.numeric(as.character(dat.n$age))
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
ggboxplot(dat.n, x = "cluster", y = "age",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text = element_text(color="black",size = 16,face="bold"),axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Age (years)")+ theme(legend.position = "none") 

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
dat.n2=dat.n[dat.n$age!="N/A",]
dat.n2$age=as.numeric(as.character(dat.n2$age))
mu <- ddply(dat.n2, "cluster", summarise, grp.mean=mean(age))
dat.n2$cluster=factor(dat.n2$cluster, levels= c("UV-low","UV-high"))
colnames(dat.n2)

####Figure2C
dat.n3=dat.n2[dat.n2$cluster=="UV-low",]
p=ggplot(dat.n3, aes(x=age)) +geom_density()+theme_classic()+theme(text = element_text(size=12))+scale_color_manual(values=c("#0000FF"))+theme(axis.title=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face="bold"))+ylab("Denstiy")+ scale_x_continuous(name = "Age (years)",breaks = seq(20, 80, 20),limits=c(10, 100))+ theme(legend.position = "none")
d <- ggplot_build(p)$data[[1]]
rand1 <- runif(1,0,40)
p + geom_area(data = subset(d, x > 40), aes(x=x, y=y), fill="blue",alpha=0.5) 

dat.n3=dat.n2[dat.n2$cluster=="UV-high",]
max(dat.n3$age)
min(dat.n3$age)
p=ggplot(dat.n3, aes(x=age)) +geom_density()+theme_classic()+theme(text = element_text(size=12))+scale_color_manual(values=c("#FF0000"))+theme(axis.title=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face="bold"))+ylab("Denstiy")+ scale_x_continuous(name = "Age (years)",breaks = seq(20, 80, 20),limits=c(10, 100))+ theme(legend.position = "none")
d <- ggplot_build(p)$data[[1]]
rand1 <- runif(1,0,40)
p + geom_area(data = subset(d, x > 40), aes(x=x, y=y), fill="red",alpha=0.5) 


##Supplementary Figure S9. Signature refitting analysis in ICGC and SRA dataset 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
dat.n=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE,quote = "")
?read.table
head(dat.n)
dat.n=t(dat.n)
rownames(dat.n)
dat.n2=as.matrix(dat.n[11:30,])
dat.n2=as.matrix(dat.n2)
dat.n2=apply(dat.n2,1,as.numeric)
dat.n2=t(dat.n2)
colnames(dat.n2)=colnames(dat.n)
my_sample_col=as.data.frame(colnames(dat.n2))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n2)
dat.n=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col$cluster=dat.n$NMF_cluster
dat.n2=dat.n2[,order(dat.n$Basis2_NMF)]
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(cluster=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE)

#Cosine similarity 
dat.n=read.table("SRA_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(dat.n$Basis2_NMF),]
ggplot(dat.n2, aes(x=factor(1:length(cos.sim.1)), y=cos.sim.1)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) + theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

###ICGC COSMIC signatures signature refitting analysis visualization
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
dat.n=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE,quote = "")
dat.n=t(dat.n)
dat.n2=as.matrix(dat.n[19:38,])
dat.n2=as.matrix(dat.n2)
dat.n2=apply(dat.n2,1,as.numeric)
dat.n2=t(dat.n2)
colnames(dat.n2)=colnames(dat.n)
my_sample_col=as.data.frame(colnames(dat.n2))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n2)
dat.n=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col$cluster=dat.n$NMF_cluster
dat.n2=dat.n2[,order(dat.n$Basis2_NMF)]
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(cluster=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend=FALSE)

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
dat.n=read.table("ICGC_annotation.txt", sep="\t", row.names=1, header=TRUE,quote = "")
dat.n2=dat.n[order(dat.n$Basis2_NMF),]
ggplot(dat.n2, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) + theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

##Figure 3. Prognostic implications of mutational signatures. 
#overall survival analysis of TCGA cohort 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/04_Clinical")
dat.n = read.table("survival_excluded.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
colnames(dat.n)
dat.n2=dat.n[complete.cases(dat.n$fustat),]
median(dat.n2$futime)
surv_object <- Surv(time = dat.n2$futime, event = dat.n2$fustat)
fit1 <- survfit(surv_object ~ cluster, data = dat.n2)
#Figure3a
ggsurvplot(fit1, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
ggsurvplot(fit1, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE,risk.table=T, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))

###Univariate analysis of overall survival 
#table 1
fit.coxph <- coxph(surv_object ~  cluster, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  age, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  Mutation_cluster2, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  sex, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)
dat.n3= dat.n2[dat.n2$stage!=c("N/A"),]
surv_object2 <- Surv(time = dat.n3$futime, event = dat.n3$fustat)
fit.coxph2 <- coxph(surv_object2 ~  stage, data = dat.n3)
ggforest(fit.coxph2, data = dat.n3,fontsize	=1, refLabel="REF",noDigits	=3)

###Multivariate analysis of overall survival 
#table 1
covariate_names <- c(cluster="cluster",age="age",Mutation_type="Mutation_cluster2",Stage="stage",sex="sex")
result=analyse_multivariate(dat.n3, vars(futime, fustat),covariates = vars(cluster, age, sex, Mutation_cluster2, stage), covariate_name_dict = covariate_names)
forest_plot(result, factor_labeller = covariate_names,  endpoint_labeller = c(futime="OS"),  orderer = ~order(HR),  labels_displayed = c("endpoint", "factor"),    ggtheme = ggplot2::theme_bw(base_size = 12))

#Landmark Survival analysis: TCGA dataset
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/04_Clinical")
dat.n = read.table("survival_excluded.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
colnames(dat.n)
dat.n2=dat.n[complete.cases(dat.n$fustat),]
median(dat.n2$futime)
dat.n2=dat.n2[dat.n2$futime>=365,]
dat.n3 = dat.n2 %>%  mutate(lm_futime = futime - 365)
lm_fit <- survfit(Surv(dat.n3$lm_futime,  dat.n3$fustat) ~ cluster, data = dat.n3)
fit.coxph=coxph(Surv(dat.n2$futime,  dat.n2$fustat) ~ cluster, subset = futime >= 365,  data = dat.n2)
#Figure 1c
ggforest(fit.coxph, data = dat.n3,fontsize	=1, refLabel="REF",noDigits	=3)
ggsurvplot(lm_fit, data = dat.n3, pval = F, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
ggsurvplot(lm_fit, data = dat.n3, pval = F, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE,risk.table=T, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))

###Univariate analysis of landmark Survival at 1-year
surv_object=Surv(dat.n2$futime,  dat.n2$fustat)
fit.coxph <- coxph(surv_object ~ cluster, subset = futime >= 365,  data = dat.n2)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  age, subset = futime >= 365,  data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  Mutation_cluster2,  subset = futime >= 365, data = dat.n2)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  sex,  subset = futime >= 365, data = dat.n2)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
dat.n3= dat.n2[dat.n2$stage!=c("N/A"),]
fit.coxph2 <- coxph(Surv(time = dat.n3$futime, event = dat.n3$fustat) ~  stage, data = dat.n3,  subset = futime >= 365)
ggforest(fit.coxph2, data = dat.n3,fontsize	=1, refLabel="REF",noDigits	=3)

###Multivariate analysis of landmark Survival at 1-year
covariate_names <- c(cluster="cluster",age="age",Mutation_type="Mutation_cluster2",Stage="stage",sex="sex")
result=analyse_multivariate(dat.n3, vars(futime, fustat),covariates = vars(cluster, age, Mutation_cluster2, stage), covariate_name_dict = covariate_names)
forest_plot(result, factor_labeller = covariate_names,  endpoint_labeller = c(futime="OS"),  orderer = ~order(HR),  labels_displayed = c("endpoint", "factor"),    ggtheme = ggplot2::theme_bw(base_size = 12))

#overall survival analysis of ICGC cohort 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/03_Validation")
dat.n = read.table("ICGC_annotation.txt", header=T, sep="\t", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
dat.n2=dat.n[complete.cases(dat.n$fustat),]
median(dat.n2$futime)
surv_object <- Surv(time = dat.n2$futime, event = dat.n2$fustat)
fit1 <- survfit(surv_object ~ NMF_cluster, data = dat.n2)
table(dat.n2$NMF_cluster)
#Figure 3b
ggsurvplot(fit1, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE, risk.table	=F,ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
ggsurvplot(fit1, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE, risk.table	=T,ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
###Univariate analysis of overall survival 
fit.coxph <- coxph(surv_object ~  NMF_cluster, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)

#Landmark Survival analysis: ICGC dataset
dat.n = read.table("ICGC_annotation.txt", header=T, sep="\t", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
dat.n2=dat.n[complete.cases(dat.n$fustat),]
median(dat.n2$futime)
dat.n2=dat.n2[dat.n2$futime>=365,]
dat.n3 = dat.n2 %>%  mutate(lm_futime = futime - 365)
lm_fit <- survfit(Surv(dat.n3$lm_futime,  dat.n3$fustat) ~ NMF_cluster, data = dat.n3)
fit.coxph=coxph(Surv(dat.n2$futime,  dat.n2$fustat) ~ NMF_cluster, subset = futime >= 365,  data = dat.n2) 
ggforest(fit.coxph, data = dat.n3,fontsize	=1, refLabel="REF",noDigits	=3)
#Figure 3d
ggsurvplot(lm_fit, data = dat.n3, pval = F, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
ggsurvplot(lm_fit, data = dat.n3, pval = F, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#FF0000","#0000FF"),conf.int = TRUE, risk.table	=T,ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 14,face="bold"),axis.text.y = element_text(size = 14, face="bold"), axis.title.x = element_text(size = 17,face="bold"),axis.title.y = element_text(size=17, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))

#Supplementary Figure S10. Clinical significance of mutational signatures in TCGA cohort 
#Figure S10a
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/04_Clinical")
dat.n = read.table("survival_excluded.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
dat.n2= dat.n[dat.n$stage!=c("N/A"),]
dat.n3= dat.n2[dat.n2$stage==c("Stage I/II"),]
dat.n3= dat.n2[dat.n2$stage==c("Stage III/IV"),]
surv_object6 <- Surv(time = dat.n3$futime, event = dat.n3$fustat)
fit6 <- survfit(surv_object6 ~ cluster, data = dat.n3)
ggsurvplot(fit6, data = dat.n3, pval = FALSE, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14),legend.position="none"))

fit.coxph <- coxph(surv_object6 ~  cluster, data = dat.n3)
ggforest(fit.coxph, data = dat.n3,fontsize	=1, refLabel="REF",noDigits	=3)

dat.n = read.table("survival_excluded.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
dat.n2=dat.n[dat.n$age_binary==c("Age > 40"),]
dat.n2=dat.n[dat.n$age_binary==c("Age < 40"),]
dat.n2=dat.n[dat.n$Mutation_cluster2==c("BRAF"),]
dat.n2=dat.n[dat.n$Mutation_cluster2==c("Non-BRAF"),]
dat.n2=dat.n[dat.n$sex==c("female"),]
dat.n2=dat.n[dat.n$sex==c("male"),]
surv_object6 <- Surv(time = dat.n2$futime, event = dat.n2$fustat)
fit6 <- survfit(surv_object6 ~ cluster, data = dat.n2)
ggsurvplot(fit6, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("UV-high", "UV-low"),xscale=365, break.x.by=730,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
fit.coxph <- coxph(surv_object6 ~  cluster, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)

#FigureS10b
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
dat.n2=dat.n[dat.n$Stage_binary!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$Stage_binary))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_text(size = 14,face="bold"),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")
chisq.test(table(dat.n$Stage_binary, dat.n$cluster,exclude="N/A"))

dat.n2=dat.n[dat.n$Ulcer!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$Ulcer))
dat.n2$Var2=factor(dat.n2$Var2,levels=c("No ulcer","Ulcer"))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title = element_text(size = 14,face="bold"),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")
chisq.test(table(dat.n$Ulcer, dat.n$cluster,exclude="N/A"))

dat.n2=dat.n[dat.n$clark_binary!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$clark_binary))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title = element_text(size = 14,face="bold"),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")
chisq.test(table(dat.n$clark_binary, dat.n$cluster,exclude="N/A"))

##Figure 4. Genomic and functional characterization. 
#Mutation analysis 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
MAF_annot=read.table("maf_annotation.csv", sep=",", header=TRUE, row.names =1, stringsAsFactors	=FALSE)
MAF_annot=as.data.frame(MAF_annot)
laml = read.maf(maf = "somatic.maf",clinicalData=MAF_annot)
?read.maf
getSampleSummary(laml)
summary = getSampleSummary(laml)
write.csv(x=summary, file="maf_summary.csv", row.names=TRUE)

#Mutation analysis result visualization
MAF_annot=read.table("maf_annotation.csv", sep=",", header=TRUE, row.names =1,stringsAsFactors	=FALSE)
MAF_annot2=MAF_annot[1:466,]
laml2=subsetMaf(laml, tsb = MAF_annot2$Tumor_Sample_Barcode, query="Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site','Frame_Shift_Del', 'Frame_Shift_Ins')")
annot_order= as.character(MAF_annot2$order)
fabcolors = c("#FF0000","#0000FF")
names(fabcolors) = c("UV-low", "UV-high")
fabcolors = list(FAB_classification = fabcolors)

#Figure 4a
oncoplot(subsetMaf(laml2, tsb = MAF_annot2$Tumor_Sample_Barcode, query="Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site')"), genes = c("BRAF", "NRAS", "HRAS", "KRAS", "NF1", "CDKN2A", "TP53", "PPP6C", "ARID2", "PTEN", "IDH1", "MAP2K1", "DDX3X", "RAC1", "RB1","KIT"), keepGeneOrder=TRUE,clinicalFeatures=c("order","cluster"),drawRowBar=FALSE, drawColBar=FALSE,removeNonMutated	=FALSE, sortByAnnotation=TRUE,groupAnnotationBySize=FALSE,fill=TRUE,annotationOrder=annot_order, annotationColor=fabcolors,bgCol="#F5F5F5")

#Figure 4b
df2 <- data.frame(supp=rep(c("BRAF", "RAS", "NF1", "Triple WT","KIT"), 2),
                  dose=rep(c("UV-low", "UV"),each=5),
                  len=c(22.7,13.3,5.3,58.7,12.0,48.2,29.5,9.5,12.8,2.6))
##UV-low 75:BRAF 17 RAS 10 NF1 4 TWT 44 KIT 9
##UV-high 390: BRAF 188 RAS 115 NF1 37 TWT 50 KIT10
df2$dose=factor(df2$dose,levels=c("UV-low","UV-high"))
df2$supp=factor(df2$supp,levels=c("BRAF", "RAS", "NF1", "Triple WT","KIT"))

ggplot(data=df2, aes(x=supp, y=len, fill=dose)) + geom_bar(stat="identity", position=position_dodge())+  geom_text(aes(label=len), vjust=1.2, color="white", position = position_dodge(0.9), size=5)+  scale_fill_manual(values=c("#0000FF","#FF0000"))+  theme_classic() +theme(axis.text = element_text( size = 14,face="bold",color="black"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion (%)")+ theme(legend.position = "none") 

##copy number alteration analysis
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/05_DNA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))

#Figure 4c
ggboxplot(dat.n, x = "cluster", y = "NEXUS_Genome",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 16,face="bold",color="black"),axis.text.y = element_text(size = 16), axis.title.y = element_text( size=18, face="bold")) + theme(legend.position = "none") + ylab("Genome with CNA (%)")
wilcox.test(NEXUS_Genome ~ cluster, data = dat.n)

#Figure 4d
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/05_DNA")
ak=read.table("UV_low_seg.csv", header=TRUE, sep=",")
ak=as.data.frame(ak)
plotFreq(ak,thres.gain=0.2,cex.lab=1.2, cex.axis=1.0)

ak=read.table("UV_high_seg.csv", header=TRUE, sep=",")
ak=as.data.frame(ak)
plotFreq(ak,thres.gain=0.2,cex.lab=1.2, cex.axis=1.0)

###DEG analysis
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/06_RNA")
x=read.table("filtered_RNA_raw.txt", sep="\t", row.names=1, header=TRUE)
group <- read.table("RNA_annotation.txt", sep="\t", row.names=1, header=FALSE)
group$V2=factor(group$V2,levels=c("UV-low", "UV-high"))
y <- DGEList(counts=x, group=group$V2)
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 117.25)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=ncol(fit$design), contrast=c(1,-1))
DEG.edgeR <- topTags(qlf, n=Inf, sort.by = "logFC", adjust.method = "BH")$table
DEG.edgeR$gene <- rownames(DEG.edgeR)
DEG.edgeR_A <- filter(DEG.edgeR, abs(logFC) >= 2 & FDR < 0.01)
write.csv(x=DEG.edgeR_A, file="DEG.csv", row.names=TRUE)

###Volcano plot 
volcanoData <- cbind(DEG.edgeR$logFC, -log10(DEG.edgeR$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
rownames(volcanoData)=rownames(DEG.edgeR)
DEGs = DEG.edgeR$FDR < 0.01 & abs(DEG.edgeR$logFC) > 2
point.col <- ifelse(DEGs, "red", "black")
mutateddf <- mutate(DEG.edgeR, significant=ifelse(DEG.edgeR$FDR < 0.01 & abs(DEG.edgeR$logFC) > 2, "Sig", "Not Sig"))
input <- cbind(mutateddf) 
input$order=abs(input$logFC)*-log10(input$FDR)
input[input$significant=="Not Sig",]$order=0
input=input[order(input$order,decreasing=TRUE),]
#Figure 4e
ggplot(input, aes(logFC, -log10(FDR))) + geom_point(aes(col=significant)) +  scale_color_manual(values=c("black", "red")) + theme(text = element_text(size=12))+geom_label_repel(data=head(input,15), aes(label=gene))+ theme_classic() +theme(axis.text.x = element_text(size=14,face="bold"),axis.text.y = element_text(size=14, face="bold"), axis.title = element_text(size=16, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14))+ theme(legend.position = "none")

###Tumor immunity markers
#Figure 4f
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
ggboxplot(dat.n, x = "cluster", y = "purity",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none") +ylab("Tumor purity")+ theme(legend.position = "none")
wilcox.test(purity ~ cluster, data = dat.n)

ggboxplot(dat.n, x = "cluster", y = "Leukocyte",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA) + theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+ylab("Leukocyte fraction")+ theme(legend.position = "none")
wilcox.test(Leukocyte ~ cluster, data = dat.n)

dat.n2=dat.n[is.na(dat.n$CYT_score_RSEM)!=TRUE,]
dat.n2$CYT_score_RSEM=as.numeric(as.character(dat.n2$CYT_score_RSEM))
ggboxplot(dat.n2, x = "cluster", y = "CYT_score_RSEM",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA) + theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+ylab("CYT score (RSEM)")+ theme(legend.position = "none")+scale_y_continuous(trans='log10')
wilcox.test(CYT_score_RSEM ~ cluster, data = dat.n2)

dat.n2=dat.n[is.na(dat.n$TIL_percentage)!=TRUE,]
dat.n2$TIL_percentage=as.numeric(as.character(dat.n2$TIL_percentage))
ggboxplot(dat.n2, x = "cluster", y = "TIL_percentage",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+ylab("Proportion of tumor infiltrating leukocytes")+theme(legend.position = "none")
wilcox.test(TIL_percentage ~ cluster, data = dat.n2)

dat.n2=dat.n[is.na(dat.n$TCR)!=TRUE,]
dat.n2$TCR=as.numeric(as.character(dat.n2$TCR))
ggboxplot(dat.n2, x = "cluster", y = "TCR",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+  theme(legend.position = "none") +ylab("TCR richness")+ theme(legend.position = "none")+scale_y_continuous(trans='log10')
wilcox.test(TCR ~ cluster, data = dat.n2)

###Cibersortx differential test 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/06_RNA")
dat.n=read.table("cibersortx.csv", sep=",", quote="",row.names=1, header=TRUE)
dat.n2=as.matrix(t(dat.n)[1:22,])
p_values <- vector("list", nrow(dat.n2))
for(i in seq_along(1: nrow(dat.n2))){
  p_values[i] = wilcox.test(dat.n2[i,1:75],dat.n2[i,76:465], paired = FALSE, alternative = "two.sided")$p.value
}
p_values = data.frame(p_values = sapply(p_values, c))
rownames(p_values)=rownames(dat.n2)
write.table(x=p_values, file="cibersortx_pvalue.txt",quote = FALSE, sep="\t",row.names = T,col.names = T)


#Supplementary Figure S11. In silico immunoprofiling analysis from RNA expression profiles using CIBERSORTx tool. 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/06_RNA")
dat.n=read.table("cibersortx.csv", sep=",", quote="",row.names=1, header=TRUE)
dat.n$cluster=rep(c("UV-low", "UV-high"), c(75,390))
dat.n=as.data.frame(dat.n)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))

ggboxplot(dat.n, x = "cluster", y = "Absolute",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")+ylab("Total immune cell fraction")
wilcox.test(Absolute ~ cluster, data = dat.n)

ggboxplot(dat.n, x = "cluster", y = "Macrophages_M1",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
wilcox.test(Macrophages_M1 ~ cluster, data = dat.n)

ggboxplot(dat.n, x = "cluster", y = "T_cells_CD4_memory_activated",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
wilcox.test(T_cells_CD4_memory_activated ~ cluster, data = dat.n)

ggboxplot(dat.n, x = "cluster", y = "T_cells_regulatory",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
wilcox.test(T_cells_regulatory ~ cluster, data = dat.n)

ggboxplot(dat.n, x = "cluster", y = "T_cells_CD8",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
wilcox.test(T_cells_CD8 ~ cluster, data = dat.n)


##Supplementary Figure S12. Comparison of the signature-based clustering with known subtypes from TCGA study. 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))

dat.n2=dat.n[dat.n$tcga_uv!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$tcga_uv))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_blank(),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")

dat.n2=dat.n[dat.n$type!="BRAF/RAS",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$type))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_blank(),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")

dat.n2=dat.n[dat.n$tcga_RNA!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$tcga_RNA))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_blank(),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")


dat.n2=dat.n[dat.n$tcga_mir!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$tcga_mir))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_blank(),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")

dat.n2=dat.n[dat.n$tcga_methyl!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$tcga_methyl))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_blank(),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")

dat.n2=dat.n[dat.n$tcga_protein!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$tcga_protein))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_blank(),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")


setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Sample_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))

table(dat.n$tcga_uv, dat.n$cluster,exclude="N/A")
chisq.test(matrix(c(37,4,2,217),nrow = 2, ncol = 2))

table(dat.n$type, dat.n$cluster,exclude="BRAF/RAS")
chisq.test(matrix(c(17,188,58,202),nrow = 2, ncol = 2))
chisq.test(matrix(c(4,37,71,353),nrow = 2, ncol = 2))
chisq.test(matrix(c(10,115,65,275),nrow = 2, ncol = 2))
chisq.test(matrix(c(44,50,31,340),nrow = 2, ncol = 2))

table(dat.n$tcga_RNA, dat.n$cluster,exclude="N/A")
chisq.test(matrix(c(13,129,25,101),nrow = 2, ncol = 2))
chisq.test(matrix(c(21,58,17,172),nrow = 2, ncol = 2))
chisq.test(matrix(c(4,43,34,187),nrow = 2, ncol = 2))


table(dat.n$tcga_mir, dat.n$cluster,exclude="N/A")
chisq.test(matrix(c(8,62,31,162),nrow = 2, ncol = 2))
chisq.test(matrix(c(12,60,27,164),nrow = 2, ncol = 2))
chisq.test(matrix(c(10,59,29,165),nrow = 2, ncol = 2))
chisq.test(matrix(c(9,43,60,181),nrow = 2, ncol = 2))

table(dat.n$tcga_methyl, dat.n$cluster,exclude="N/A")
chisq.test(matrix(c(10,54,30,177),nrow = 2, ncol = 2))
chisq.test(matrix(c(17,62,23,169),nrow = 2, ncol = 2))
chisq.test(matrix(c(4,69,36,162),nrow = 2, ncol = 2))
chisq.test(matrix(c(9,46,31,185),nrow = 2, ncol = 2))

table(dat.n$tcga_protein, dat.n$cluster,exclude="N/A")
chisq.test(matrix(c(4,10,14,97),nrow = 2, ncol = 2))
chisq.test(matrix(c(5,45,13,101),nrow = 2, ncol = 2))
chisq.test(matrix(c(5,42,13,104),nrow = 2, ncol = 2))
chisq.test(matrix(c(4,10,14,136),nrow = 2, ncol = 2))

###panel sequencing data preparation
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel")
sig=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
sig=as.data.frame(sig)
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)
head(sig)
sample=rownames(sig)
#Deconstrucsig of MSK panel sequencing data 
for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = list1)
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/COSMIC_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/COSMIC_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig2 = read.table("deconstruct.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(signatures.exome.cosmic.v3.may2019)
for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =sig2, contexts.needed = FALSE, associated	=c("SigA", "SigB"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/basis_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/basis_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}


###PANEL simulation data preparation
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
dat.n=read.table("simulate_deconstruc.txt", sep="\t",  header=TRUE)
sim.sig <- mut.to.sigs.input(mut.ref = dat.n, sample.id = "Sample", chr = "chr",pos = "pos",  ref = "ref", alt = "alt",bsg = BSgenome.Hsapiens.UCSC.hg38)
sample=dat.n$Sample
for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sim.sig, signatures.ref =signatures.exome.cosmic.v3.may2019, contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"),  sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/SIM_COSMIC_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/SIM_COSMIC_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/SIM_COSMIC_TRI")
  write.table(x=tumor, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig2 = read.table("deconstruct.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(signatures.exome.cosmic.v3.may2019)
for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sim.sig, signatures.ref =sig2, contexts.needed = TRUE, associated	=c("SigA", "SigB"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/SIM_basis_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel/SIM_basis_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

##Figure 5. Independent panel sequencing cohort analysis and in silico panel simulation of TCGA cohortVisualization  of MSK cohort (panel sequencing cohort)
#Figure 5a. Visualization  of MSK cohort (panel seqeuncing cohort)
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig2 = read.table("deconstruct.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(signatures.exome.cosmic.v3.may2019)
rownames(sig2)
dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = "MSK_UV_low")
pdf('MSK_UV_low.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()

dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = "MSK_UV_high")
pdf('MSK_UV_high.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()


#Figure 5b. Visualization  of TCGA cohort (in silico panel simulation)
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/02_Deconstruct")
sig2 = read.table("deconstruct.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(signatures.exome.cosmic.v3.may2019)
dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id =  "SIM_UV_low" )

pdf('SIM_UV_low.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()

dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id =  "SIM_UV_high" )
pdf('SIM_UV_high.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()

#Figure 5c, upper 
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel")
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(dat.n$SigB),]
ggplot(dat.n2, aes(x=factor(1:length(basis_cos)), y=as.numeric(as.character(basis_cos)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n=as.data.frame(dat.n[,6:8])
dat.n=t(dat.n)
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Cluster"
dat.n=dat.n[,order(sample$SigB)]
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Cluster=Var2)
##1200/600
pheatmap(dat.n,show_rownames=TRUE,show_colnames=FALSE,cluster_cols=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col =my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)

#Figure 5c, lower
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(dat.n$SigB),]
ggplot(dat.n2, aes(x=factor(1:length(COSMIC_cos)), y=as.numeric(as.character(COSMIC_cos)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(as.data.frame(dat.n[,12:31]))
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Cluster"
dat.n=dat.n[,order(sample$SigB)]
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list( Cluster=Var2)
pheatmap(dat.n,show_rownames=TRUE,show_colnames=FALSE,cluster_cols=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col =my_sample_col,color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,fontsize = 12,legend=F,annotation_legend = F)

#Figure 5d
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
ggboxplot(dat.n, x = "cluster", y = "mutation",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text = element_text( size = 14,face="bold",color="black"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
wilcox.test(mutation ~ cluster, data = dat.n)

ggboxplot(dat.n, x = "cluster", y = "Sum_SBS7",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text = element_text( size = 14,face="bold",color="black"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Sum of SBS7")+ theme(legend.position = "none") 
wilcox.test(Sum_SBS7 ~ cluster, data = dat.n)



###Figure S13. Comparison of mutational context of WES and in silico panel simulation in TCGA cohort. 
#Figure S13a
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel")
dat.n3=read.table("simulate_annotation.txt", sep="\t", row.names=1, header=TRUE,stringsAsFactors=FALSE)
dat.n3=dat.n3[order(dat.n3$orig.order),]
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA")
dat.n=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
rownames(dat.n)
dat.n=dat.n[rownames(dat.n3),]
dat.n=t(dat.n)
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="Cluster"
rownames(my_sample_col)=my_sample_col$Cluster
my_sample_col$Cluster=rep(c("UV-low","UV-high"), c(13,370))
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
dat.n = dat.n[,order(dat.n2$SigB)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE)

dat.n3=t(dat.n)
sum=cbind(rowSums(dat.n3[,1:16]),rowSums(dat.n3[,17:32]),rowSums(dat.n3[,33:48]),rowSums(dat.n3[,49:64]),rowSums(dat.n3[,65:80]),rowSums(dat.n3[,81:96]))
sum=as.data.frame(sum)
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% tidyr::gather(variable, value, -ind)
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")

#Figure S13b
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel")
dat.n=read.table("simulate_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
dat.n=dat.n[order(dat.n$orig.order),]
median(dat.n$mut.retention)
ggplot(dat.n, aes(x=factor(1:length(mut.retention)), y=as.numeric(as.character(mut.retention)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel")
dat.n=read.table("simulate_tri.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=read.table("simulate_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
dat.n2=dat.n[order(dat.n2$orig.order),]
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/")
dat.n=read.table("Heatmap.txt", sep="\t", row.names=1, header=TRUE)
dat.n3=dat.n[rownames(dat.n2),]
dat.n2=as.matrix(dat.n2)
dat.n3=as.matrix(dat.n3)
cos=cos.sim(dat.n2,dat.n3)
cos.matrix=as.data.frame(diag(cos))
colnames(cos.matrix)="cos.sim"
median(cos.matrix$cos.sim)
ggplot(cos.matrix, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

dat.n2=as.data.frame(t(dat.n2))
dat.n3=as.data.frame(t(dat.n3))
my_sample_col=as.data.frame(colnames(dat.n2))
colnames(my_sample_col)="Cluster"
rownames(my_sample_col)=colnames(dat.n2)
my_sample_col$Cluster=rep(c("UV-low","UV-high"), c(13,370))
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n2)
colnames(my_gene_col)="Type"
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low","UV-high")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend=TRUE)

dat.n2=t(dat.n2)
sum=as.data.frame(cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96])))
sum=as.data.frame(sum)
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% tidyr::gather(variable, value, -ind)
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")


###Figure S14. Panel sequencing cohort (MSK cohort) classified by dominant signatures.
setwd("/data/Delta_data4/kysbbubbu/Frontier_TCGA/07_Panel")
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
dat.n=dat.n[order(dat.n$SigB),]
ggplot(dat.n, aes(x=factor(1:length(mutation)), y=mutation)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

dat.n=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Cluster"
dat.n=dat.n[,order(sample$SigB)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("UV-low", "UV-high")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n,show_rownames=FALSE,show_colnames=FALSE,cluster_cols=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col =my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)

dat.n=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
dat.n2=dat.n[order(sample$SigB),]
sum=as.data.frame(cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96])))
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% tidyr::gather(variable, value, -ind)
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")