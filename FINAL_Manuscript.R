##Library Load
library(devtools)
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

###Estimate the rank of NMF decomposition 
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
NMF=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
NMF = as.matrix(NMF)
NMF = t(NMF)
estim.r <- nmf(NMF, 2:8, nrun=1000, seed=123456, .options="p24")
?nmf
#1000/700
plot(estim.r)
summary(estim.r)

#1200/600
consensusmap(estim.r, labCol=NA, labRow=NA,tracks = c("consensus:"))


###Figure S1
cophenetic=as.data.frame(estim.r$measures$cophenetic)
colnames(cophenetic)="cophenetic"
cophenetic$rank=estim.r$measures$rank
cophenetic$rss=estim.r$measures$rss
cophenetic=cophenetic[1:5,]
#500/400
ggplot(cophenetic, aes(x=factor(rank), y=cophenetic,group=1)) +geom_line()+geom_point(color="black") + theme_classic()+theme(axis.title=element_text(size=16,face="bold"),axis.text=element_text(size=14),plot.background=element_blank(), axis.line = element_line(colour="black"))+labs(x = "NMF Rank", y = "Cophenetic coefficient")


###Figure S1
cophenetic=as.data.frame(estim.r$measures$silhouette.coef)
colnames(cophenetic)="cophenetic"
cophenetic$rank=estim.r$measures$rank
cophenetic=cophenetic[1:5,]
#500/400
ggplot(cophenetic, aes(x=factor(rank), y=cophenetic,group=1)) +geom_line()+geom_point(color="black") + theme_classic()+theme(axis.title=element_text(size=16,face="bold"),axis.text=element_text(size=14),plot.background=element_blank(), axis.line = element_line(colour="black"))+labs(x = "NMF Rank", y = "Average silhouette width")


#NMF decomposition using the rank (2)
NMF=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
NMF = t(NMF)
res <- nmf(NMF, 2, nrun=1000, seed=123456, .options="p24")
res
consensusmap(res, labCol=NA, labRow=NA,tracks = c("consensus:"))
fit(res)
res
basismap(res,Rowv=NA, Colv=NA)
coefmap(res,Rowv=NA, Colv=TRUE,distfun="euclidean",hclustfun="complete")
basis=as.matrix(res@fit@W)
basis=t(t(basis)/colSums(basis))
pheatmap(basis,cluster_cols=FALSE, show_rownames=FALSE, labels_col=c("SigA","SigB"), angle_col=0, color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),border_color=NA,  fontsize_col=14, cellwidth=130)
##write.csv(x=basis, file="basis.csv", row.names=TRUE)
coefmap=as.matrix(res@fit@H)
coefmap=t(t(coefmap)/colSums(coefmap))
pheatmap(coefmap,cluster_cols=TRUE, show_rownames=TRUE,labels_row	=c("SigA","SigB"),show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),clustering_method ="ward.D2")

#K-means clustering 
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=cbind(dat.n$SigA,dat.n$SigB)
dat.n2=as.data.frame(dat.n2)
dat.n2$cluster=kmeans(dat.n2[,1:2],2)$cluster
ggplot2.scatterplot(data=dat.n2, xName='V1',yName='V2', groupName='cluster') + labs(x = "SigA", y="SigB")+ theme_classic()+theme(text = element_text(size=14))+theme(axis.text.x = element_text( size = 12,face="bold"),axis.text.y = element_text(size = 12,face="bold"), axis.title.x = element_text( size=12, face="bold"),axis.title.y = element_text( size=12, face="bold"))+  scale_colour_manual(name="cluster", values= c("red", "blue"))+ theme(legend.position = "none") 
dat.n2$cluster=as.factor(dat.n2$cluster)
current.cluster.ids <- c("1","2")
new.cluster.ids <- c("UV","Non-UV")
dat.n2$cluster  <- plyr::mapvalues(x = dat.n2$cluster, from = current.cluster.ids, to = new.cluster.ids)
dat.n2$cluster=factor(dat.n2$cluster, levels= c("UV-low","UV-high"))
NMF_cluster=dat.n2$cluster
summary(NMF_cluster)
#write.csv(x=NMF_cluster, file="cluster.csv", row.names=TRUE)


dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
colnames(dat.n)
dat.n2=dat.n[order(dat.n$SigB),]
ggplot(dat.n2, aes(x=factor(1:length(Mutation)), y=Mutation)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))
?

###Figure 1C
coefmap2 = coefmap[,order(dat.n$SigB)]
rownames(coefmap2)=c("SigA","SigB")
my_gene_col <- data.frame(sample = rep(c("SigA", "SigB"), c(1,1)))
rownames(my_gene_col)=rownames(coefmap2)
colnames(my_gene_col)="Basis"
my_gene_col$Basis=factor(my_gene_col$Basis, levels= c("SigA", "SigB"))
pheatmap(coefmap2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_names_row=FALSE,annotation_row =my_gene_col,annotation_names_col=FALSE,annotation_colors=ann_colors)

dat.n=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n2=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
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
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Type = Var1, Cluster=Var2)
##900/600
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 11, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 11, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = TRUE,legend	=T,cellwidth	=0.6)



###Figure S2
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, sample.id = "Sample", chr = "chr",pos = "pos",ref = "ref",alt = "alt")
sig = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig=as.data.frame(t(sig))
colnames(sig)=colnames(sigs.input)
rownames(sig)

##cos 0.954
SBS7_low = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'SBS7_low')
?whichSignatures
tumor=SBS7_low$tumor
deconstruct=SBS7_low$product
cos=cos.sim(tumor,deconstruct)
cos

##cos 0.999
SBS7_high = whichSignatures(tumor.ref = sig, signatures.ref= signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, sample.id = 'SBS7_high', associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
tumor=SBS7_high$tumor
deconstruct=SBS7_high$product
cos=cos.sim(tumor,deconstruct)
cos
cos.sim()

##cos 0.996
Total = whichSignatures(tumor.ref = sig, signatures.ref= signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, sample.id = 'Total', associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
tumor=Total$tumor
deconstruct=Total$product
cos=cos.sim(tumor,deconstruct)
cos
cos.sim(Total$tumor, SBS7_high$tumor)

###Figure 1A 
##cos 0.952
SigA = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, sample.id = 'SigA',associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
tumor=SigA$tumor
deconstruct=SigA$product
cos=cos.sim(tumor,deconstruct)
cos
##cos 0.999
SigB= whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019, contexts.needed = FALSE, sample.id = 'SigB',associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"))
tumor=SigB$tumor
deconstruct=SigB$product
cos=cos.sim(tumor,deconstruct)
cos
rownames(sig)
cos.sim(SigA$tumor,SBS7_low$tumor)
cos.sim(SigB$tumor,SBS7_high$tumor)

pdf('SBS7_high.pdf', width = 10, height = 10)
chart<-plotSignatures(SBS7_high)
dev.off()


pdf('SBS7_low_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SBS7_low)
dev.off()

pdf('Total_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(Total)
dev.off()

pdf('SigA_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SigA)
dev.off()


pdf('SigB_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(SigB)
dev.off()


#comparison between COSMIC SBS signatures and NMF-decomposed basis

sig=as.data.frame(t(basis))
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)
COSMIC.v3=signatures.exome.cosmic.v3.may2019

#Define cosine similarity
cos.sim=function(ma, mb){
  mat=tcrossprod(ma, mb)
  t1=sqrt(apply(ma, 1, crossprod))
  t2=sqrt(apply(mb, 1, crossprod))
  mat / outer(t1,t2)
}

#SigA
SigA_cos=as.matrix(sig[1,])
deconstruct=as.matrix(COSMIC.v3)
cos=cos.sim(SigA_cos,deconstruct)
cos=cos[,order(cos,decreasing=TRUE)]
cos=round(cos,2)
cos

#SigB
SigB_cos=as.matrix(sig[2,])
deconstruct=as.matrix(COSMIC.v3)
cos2=cos.sim(SigB_cos,deconstruct)
cos2=cos2[,order(cos2,decreasing=TRUE)]
cos2=round(cos2,2)
cos2




###Figure 1A : piechart 
barplot = as.data.frame(rep(c("SBS1", "SBS3", "SBS5", "SBS7a", "SBS7b","SBS38","Unkown"),2))
colnames(barplot)="SBS"
barplot$proportion=c("0.205","0.321","0.222","0","0","0.089","0.163","0","0","0.076","0.512","0.405","0","0.007")
barplot$basis=rep(c("SigA", "SigB"),each=7)
barplot$proportion=as.numeric(as.character(barplot$proportion))
library(scales)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
barplot$SBS=factor(barplot$SBS, levels= c("SBS1","SBS3","SBS5","SBS7a","SBS7b","SBS38","Unkown"))
ggplot(data = barplot, aes(x = "", y = proportion, fill = SBS )) +   geom_bar(stat = "identity", position = position_fill()) +  blank_theme + coord_polar(theta = "y") +  facet_wrap(~ basis)  +    theme(axis.text=element_blank(),legend.position='top',legend.title=element_blank(),legend.text=element_text(size=12))

###Figure S3
#Mutational context from Mutect2 data 
dat.n=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n2=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n = dat.n[,order(dat.n2$Heatmap_order)]
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
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend	=FALSE)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend=TRUE)

#Mutational context from Varscan2 data
dat.n=read.table("varscan_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n2=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n = dat.n[,order(dat.n2$Heatmap_order)]
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
names(Var2) = c("Non-UV", "UV")
ann_colors = list(type = Var1, cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend=TRUE)

#Cosine similarity between Mutect2 and Varscan2  
dat.n=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n2=read.table("varscan_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n=as.matrix(dat.n)
dat.n2=as.matrix(dat.n2)
cos=cos.sim(dat.n,dat.n2)
cos.matrix=diag(cos)
cos.matrix=as.data.frame(cos.matrix)
min(cos.matrix)
cos.matrix2=as.data.frame(cos.matrix[colnames(coefmap[,pheatmap$tree_col[["order"]]]),])
colnames(cos.matrix2)="cos.sim"
?pheatmap
ggplot(cos.matrix2, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 
?scale_y_continuous

###Figure S4
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
my_comparisons <- list( c("UV-low","UV-high"))
ggboxplot(dat.n, x = "cluster", y = "Read_depth.1",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text= element_text(size = 15,face="bold",color="black"), axis.title.y = element_text(size=17, face="bold"))+theme(legend.title =element_text(size = 16,face="bold"),legend.text = element_text(size=16), legend.position="right")+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +ylab("Read depth (ENSEMBL exonic region)") +labs(colour="Cluster")+ theme(legend.position = "none")
+stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

###Figure S4
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
dat.n2 <- as.data.frame(table(dat.n$cluster,dat.n$Neoadjuvant.CTx.Hx))
dat.n2$Var2=factor(dat.n2$Var2,levels=c("No","Yes"))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text= element_text(size = 15,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title =element_text(size = 16,face="bold"),legend.text = element_text(size=16), legend.position="right")+labs(fill="Neoadjuvant CTx")+NoLegend()

###Figure S4
dat.n2 <- as.data.frame(table(dat.n$cluster,dat.n$Primary.tumor.diagnosed))
dat.n2$Var2=factor(dat.n2$Var2,levels=c("No","Yes"))
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text= element_text(size = 15,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title =element_text(size = 16,face="bold"),legend.text = element_text(size=16), legend.position="right")+labs(fill="Primary tumor")+NoLegend()


###NMF decomposition of SRA dataset
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dbgap_NMF=read.table("dbGaP_heatmap2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_NMF = as.matrix(dbgap_NMF)
dbgap_NMF = t(dbgap_NMF)
dbgap_estim.r <- nmf(NMF2, 2:8, nrun=1000, seed=123456, .options="p24")
dbgap_NMF=read.table("dbGaP_heatmap2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_NMF = as.matrix(dbgap_NMF)
dbgap_NMF = t(dbgap_NMF)
dbgap_res <- nmf(dbgap_NMF, 2, nrun=1000, seed=123456, .options="p24")
basis=as.matrix(dbgap_res@fit@W)
basis=t(t(basis)/colSums(basis))
pheatmap(basis,cluster_cols=FALSE, show_rownames=FALSE, labels_col=c("Basis1","Basis2"), angle_col=0, color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),border_color=NA,  fontsize_col=14, cellwidth=130)
#write.csv(x=basis, file="dbgap_basis.csv", row.names=TRUE)
dbgap_coefmap=as.matrix(dbgap_res@fit@H)
dbgap_coefmap=t(t(dbgap_coefmap)/colSums(dbgap_coefmap))
pheatmap(dbgap_coefmap,cluster_cols=TRUE, show_rownames=TRUE,labels_row	=c("Basis1","Basis2"),show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),clustering_method ="ward.D2")
#write.csv(x=dbgap_coefmap, file="dbgap_coefmap.csv",row.names=TRUE)

###NMF decomposition of ICGC dataset
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
NMF2=read.table("ICGC_WES_tri2.txt", sep="\t", row.names=1, header=TRUE)
NMF2 = as.matrix(NMF2)
NMF2 = t(NMF2)
NMF2_estim.r <- nmf(NMF2, 2:8, nrun=1000, seed=123456, .options="p24")
NMF2=read.table("ICGC_WES_tri2.txt", sep="\t", row.names=1, header=TRUE)
NMF2 = t(NMF2)
res2 <- nmf(NMF2, 2, nrun=1000, seed=123456, .options="p24")
basis2=as.matrix(res2@fit@W)
basis2=t(t(basis2)/colSums(basis2))
#write.csv(x=basis2, file="basis2_ICGC",row.names=TRUE)
icgc_coefmap=as.matrix(res2@fit@H)
icgc_coefmap=t(t(icgc_coefmap)/colSums(icgc_coefmap))
pheatmap(icgc_coefmap,cluster_cols=TRUE, show_rownames=TRUE,labels_row	=c("Basis1","Basis2"),show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),clustering_method ="ward.D2")
#write.csv(x=icgc_coefmap, file="icgc_coefmap.csv",row.names=TRUE)

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


####ICGC basis --deconstrucSig and visulization
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




#SRA basis --TCGA basis comparison 
basis=as.matrix(res@fit@W)
basis=t(t(basis)/colSums(basis))
sig=t(basis)
colnames(sig)=colnames(signatures.exome.cosmic.v3.may2019)

basis=as.matrix(dbgap_res@fit@W)
basis=t(t(basis)/colSums(basis))
dbgap_sig=t(basis)
colnames(dbgap_sig)=colnames(signatures.exome.cosmic.v3.may2019)

dbgap_Basis1_cos=t(as.matrix(dbgap_sig[1,]))
Basis1_cos=t(as.matrix(sig[1,]))
cos.sim(dbgap_Basis1_cos,Basis1_cos)

dbgap_Basis2_cos=t(as.matrix(dbgap_sig[2,]))
Basis2_cos=t(as.matrix(sig[2,]))
cos.sim(dbgap_Basis2_cos,Basis2_cos)


####SRA basis --deconstrucSig and visulization
dbgap_Basis=as.data.frame(rbind(dbgap_Basis1_cos,dbgap_Basis2_cos))
rownames(dbgap_Basis)=c("Basis1","Basis2")
colnames(dbgap_Basis)=colnames(sigs.input)

##cos 0.998
dbgap_Basis1 = whichSignatures(tumor.ref = dbgap_Basis, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'Basis1')
cos.sim(dbgap_Basis1$tumor,dbgap_Basis1$product)


##cos 0.942
dbgap_Basis2 = whichSignatures(tumor.ref = dbgap_Basis, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = 'Basis2')
cos.sim(dbgap_Basis2$tumor,dbgap_Basis2$product)


pdf('dbgap_Basis1_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(dbgap_Basis1)
dev.off()


pdf('dbgap_Basis2_sigplot.pdf', width = 10, height = 10)
chart<-plotSignatures(dbgap_Basis2)
dev.off()


###SRA dataset coef/mutation/context visulization
###Figure S5
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_annt2=dbgap_annt[order(dbgap_annt$Basis2_NMF),]
dbgap_annt3=t(dbgap_annt2[,8:9])
colnames(dbgap_annt2)
my_sample_col=as.data.frame(dbgap_annt[,10])
rownames(my_sample_col)=rownames(dbgap_annt)
colnames(my_sample_col)="Cluster"
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Cluster=Var2)
##1000/150
pheatmap(dbgap_annt3,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend=FALSE,annotation_colors=ann_colors,legend=FALSE)

?pheatmap
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_annt2=dbgap_annt[order(dbgap_annt$Basis2_NMF),]
head(dbgap_annt2)
##800/200
ggplot(dbgap_annt2, aes(x=factor(1:length(Mutation)), y=Mutation)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
coln
dbgap_annt2=dbgap_annt[order(dbgap_annt$Basis2_NMF),]
ggplot(dbgap_annt2, aes(x=factor(1:length(Mutation)), y=Mutation)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

dat.n=read.table("dbGaP_heatmap2.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
my_gene_col <- data.frame(Type = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
head(dbgap_annt)
my_sample_col=as.data.frame(dbgap_annt[,3])
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Study"
dat.n=dat.n[,order(dbgap_annt$Basis2_NMF)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = brewer.pal(7,"Set1")
names(Var2) = c("Hugo", "Krauthammer", "Liang", "Riaz", "Roh", "Shain", "Snyder")
ann_colors = list(Type = Var1, Study=Var2)
head(dat.n)

pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col=my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors)
dbgap_heatmap=pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col=my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE)

dat.n=read.table("dbGaP_heatmap2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(dbgap_annt$Basis2_NMF),]
sum=cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
?mutate
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")



###ICGC dataset coef/mutation/context visulization
###Figure S5

icgc_annt$cluster=kmeans(icgc_annt[,7:8],2)$cluster
current.cluster.ids <- c("1","2")
new.cluster.ids <- c("UV-low","UV-high")
icgc_annt$cluster <- plyr::mapvalues(x = icgc_annt$cluster, from = current.cluster.ids, to = new.cluster.ids)
icgc_annt$cluster=factor(icgc_annt$cluster, levels= c("UV-low","UV-high"))
table(icgc_annt$cluster)

icgc_annot=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE, quote = "")
icgc_annot=icgc_annot[order(icgc_annot$Basis2_NMF),]
colnames(icgc_annot)
my_sample_col=as.data.frame(icgc_annot[,11])
icgc_annot=icgc_annot[,9:10]
dat.n=t(icgc_annot)
rownames(dat.n)=c("SigA","SigB")
colnames(my_sample_col)="Cluster"
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Cluster=Var2)
rownames(my_sample_col)=colnames(dat.n)
##1000/150
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,annotation_legend = FALSE,legend=FALSE)

dat.n=read.table("ICGC_WES_tri2.txt", sep="\t", row.names=1, header=TRUE, quote = "")
dat.n=t(dat.n)
my_gene_col <- data.frame(Type = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
dbgap_annt=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col=as.data.frame(dbgap_annt[,16])
head(dbgap_annt)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Study"
dat.n2=dat.n[,order(dbgap_annt$Basis2_NMF)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = brewer.pal(2,"Set1")
names(Var2) = c("SKCM-US", "MELA-AU")

ann_colors = list(Type = Var1, Study=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col=my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)


icgc_annt=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
icgc_annt=icgc_annt[order(icgc_annt$Basis2_NMF),]
ggplot(icgc_annt, aes(x=factor(1:length(WES_SNP)), y=WES_SNP)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

dat.n=read.table("ICGC_WES_tri2.txt", sep="\t", row.names=1, header=TRUE, quote = "")
icgc_annt=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(icgc_annt$Basis2_NMF),]

sum=cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")


###Mutation count comparison
#FigureS6
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dbgap_annt=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
head(dbgap_annt)
dbgap_annt$cluster=factor(dbgap_annt$cluster, levels=c("UV-low","UV-high"))
my_comparisons <- list( c("UV-low","UV-high") )
ggboxplot(dbgap_annt, x = "cluster", y = "Mutation",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 12,face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dbgap_annt=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_annt$NMF_cluster=factor(dbgap_annt$NMF_cluster, levels=c("UV-low","UV-high"))
ggboxplot(dbgap_annt, x = "NMF_cluster", y = "WES_mutation",color = "NMF_cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "NMF_cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text( size =  16,face="bold",color="black"),axis.text.y = element_text( size = 12,face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_annt$NMF_cluster=factor(dbgap_annt$NMF_cluster, levels=c("UV-low","UV-high"))
ggboxplot(dbgap_annt, x = "NMF_cluster", y = "Mutation",color = "NMF_cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "NMF_cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text( size =  16,face="bold",color="black"),axis.text.y = element_text( size = 12,face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

##PCA
##Figure S7
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n=as.matrix(dat.n)
icgc_annt=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n_pca <- prcomp(dat.n, center = TRUE,scale. = FALSE)
fviz_pca_ind(dat.n_pca, col.ind = factor(icgc_annt$cluster),palette = c("red",  "blue"),label="none",title="")+ theme_bw()+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ theme(legend.position = "none")

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dat.n=read.table("ICGC_WES_tri2.txt", sep="\t", row.names=1, header=TRUE, quote = "")
icgc_annt=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(icgc_annt$Basis2_NMF),]
dat.n = as.matrix(dat.n)
dat.n_pca <- prcomp(dat.n, center = TRUE,scale. = FALSE)
fviz_pca_ind(dat.n_pca,       col.ind = factor(icgc_annt$NMF_cluster),palette = c("red",  "blue"),label="none",title="")+ theme_bw()+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ theme(legend.position = "none")

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dbgap_NMF=read.table("dbGaP_heatmap2.txt", sep="\t", row.names=1, header=TRUE)
dbgap_NMF = as.matrix(dbgap_NMF)
dbgap_pca <- prcomp(dbgap_NMF, center = TRUE,scale. = FALSE)
dbgap_annt=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
fviz_pca_ind(dbgap_pca,       col.ind = factor(dbgap_annt$NMF_cluster),palette = c("red",  "blue"),label="none",title="")+ theme_bw()+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ theme(legend.position = "none")

##WGS/WES comparison
#Figure S8

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dat.n=read.table("ICGC_WES_tri2.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=read.table("ICGC_WGS_tri.txt", sep="\t", row.names=1, header=TRUE)
dat.n=dat.n[1:198,]
dat.n2=dat.n2[rownames(dat.n), ]

dat.n3=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
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

##cos.matrix2=as.data.frame(cos.matrix[colnames(coefmap[,pheatmap$tree_col[["order"]]]),])
colnames(cos.matrix)="cos.sim"
median(cos.matrix$cos.sim)
?pheatmap
ggplot(cos.matrix, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

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
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Type = Var1, Cluster=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_row =my_gene_col,annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend=FALSE,annotation_colors=ann_colors,legend=FALSE)
dat.n2=t(dat.n2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend=FALSE,annotation_colors=ann_colors,legend=FALSE)

dat.n2=t(dat.n2)
sum=cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
?mutate
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")

dat.n2=t(dat.n)
sum=cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
?mutate
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")

proportion=as.data.frame(dat.n3$WES_SNP/dat.n3$WGS_SNP*100)
colnames(proportion)="proportion"
ggplot(proportion, aes(x=factor(1:length(proportion)), y=proportion)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Proportion of SNVs")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))
median(proportion$proportion)


###TCGA cohort deconstrucsig (COSMIC signature refitting)

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, sample.id = "Sample", chr = "chr",pos = "pos",ref = "ref",alt = "alt")
sig=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
sig=as.data.frame(sig)
colnames(sig)=colnames(sigs.input)
head(sig)
rowSums(sig)
sample=read.table("list.txt", sep="\t", header=FALSE,stringsAsFactors=FALSE,colClasses	= "character",quote = "")
list=sample$V1
for (list1 in list) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/02_SIG/decon_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/02_SIG/decon_sig")
  weight=dec$weights
  write.table(x=weight, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

###TCGA cohort deconstrucsig (Basis signature refitting)

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, sample.id = "Sample", chr = "chr",pos = "pos",ref = "ref",alt = "alt")
sig2 = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(sigs.input)
sig=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
sig=as.data.frame(sig)
colnames(sig)=colnames(sigs.input)
head(sig)
head(sig2)
rowSums(sig)
sample=read.table("list.txt", sep="\t", header=FALSE,stringsAsFactors=FALSE,colClasses	= "character",quote = "")
list=sample$V1
for (list1 in list) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =sig2, contexts.needed = FALSE, associated	=c("Basis1", "Basis2"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/02_SIG/basis_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/02_SIG/basis_sig")
  weight=dec$weights
  write.table(x=weight, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}


###TCGA COSMIC signatures signature refitting analysis visualization
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/02_SIG")
dat.n=read.table("COSMIC_signature_artificial.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n=dat.n[1:20,]
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$cluster=rep(c("Non-UV", "UV"), c(75,391))
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(cluster=Var2)
##1000/500
#Figure 2A
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 11, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE,fontsize = 12)
?pheatmap


###TCGA Basis signature refitting analysis visualization
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/02_SIG")
dat.n=read.table("COSMIC_signature_artificial.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
dat.n=dat.n[26:28,]
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n)
my_sample_col$cluster=rep(c("Non-UV", "UV"), c(75,391))
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
rownames(dat.n)=c("SigA","SigB","Unknown")
ann_colors = list(cluster=Var2)
##1000/600
#Figure 2A
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 11, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors,fontsize=0,legend=F,annotation_legend=F)



####TCGA Age distribution 
####Figure2B
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$age=as.numeric(as.character(dat.n$age))
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
ggboxplot(dat.n, x = "cluster", y = "age",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text = element_text(color="black",size = 16,face="bold"),axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Age (years)")+ theme(legend.position = "none") 

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
dat.n2=dat.n[dat.n$age!="N/A",]
dat.n2$age=as.numeric(as.character(dat.n2$age))
mu <- ddply(dat.n2, "cluster", summarise, grp.mean=mean(age))
dat.n2$cluster=factor(dat.n2$cluster, levels= c("UV-low","UV-high"))
colnames(dat.n2)

####Figure2C
dat.n3=dat.n2[dat.n2$cluster=="Non-UV",]
p=ggplot(dat.n3, aes(x=age)) +geom_density()+theme_classic()+theme(text = element_text(size=12))+scale_color_manual(values=c("#0000FF"))+theme(axis.title=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face="bold"))+ylab("Denstiy")+ scale_x_continuous(name = "Age (years)",breaks = seq(20, 80, 20),limits=c(10, 100))+ theme(legend.position = "none")
d <- ggplot_build(p)$data[[1]]
rand1 <- runif(1,0,40)
p + geom_area(data = subset(d, x > 40), aes(x=x, y=y), fill="blue",alpha=0.5) 

####Figure2C
dat.n3=dat.n2[dat.n2$cluster=="UV",]
max(dat.n3$age)
min(dat.n3$age)
p=ggplot(dat.n3, aes(x=age)) +geom_density()+theme_classic()+theme(text = element_text(size=12))+scale_color_manual(values=c("#FF0000"))+theme(axis.title=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face="bold"))+ylab("Denstiy")+ scale_x_continuous(name = "Age (years)",breaks = seq(20, 80, 20),limits=c(10, 100))+ theme(legend.position = "none")
d <- ggplot_build(p)$data[[1]]
rand1 <- runif(1,0,40)
p + geom_area(data = subset(d, x > 40), aes(x=x, y=y), fill="red",alpha=0.5) 


###SRA COSMIC signatures signature refitting analysis visualization
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dat.n=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE,quote = "")
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

dat.n=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col$cluster=dat.n$NMF_cluster
dat.n2=dat.n2[,order(dat.n$Basis2_NMF)]

# Specify colors
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(cluster=Var2)
dat.n2
##1000/600
#Figure S9
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE)
?pheatmap

dat.n=read.table("dbGaP_annotation2.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(dat.n$Basis2_NMF),]
dat.n2
#800/100
ggplot(dat.n2, aes(x=factor(1:length(cos.sim.1)), y=cos.sim.1)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) + theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

###ICGC COSMIC signatures signature refitting analysis visualization
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dat.n=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE,quote = "")
dat.n=t(dat.n)
dat.n2=as.matrix(dat.n[19:38,])
dat.n2=as.matrix(dat.n2)
dat.n2=apply(dat.n2,1,as.numeric)
dat.n2=t(dat.n2)
colnames(dat.n2)=colnames(dat.n)
my_sample_col=as.data.frame(colnames(dat.n2))
colnames(my_sample_col)="cluster"
rownames(my_sample_col)=colnames(dat.n2)

dat.n=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE)
my_sample_col$cluster=dat.n$NMF_cluster
dat.n2=dat.n2[,order(dat.n$Basis2_NMF)]
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(cluster=Var2)

##1000/500
#Figure S9
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_col=my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend=FALSE)

?pheatmap

dat.n=read.table("ICGC_META_data2.txt", sep="\t", row.names=1, header=TRUE,quote = "")
dat.n2=dat.n[order(dat.n$Basis2_NMF),]
#800/100
ggplot(dat.n2, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) + theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))


#Survival analysis: ICGC dataset
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
dat.n = read.table("ICGC_META_data.txt", header=T, sep="\t", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
colnames(dat.n)
dat.n2=dat.n[complete.cases(dat.n$fustat),]
surv_object <- Surv(time = dat.n2$futime, event = dat.n2$fustat)
fit1 <- survfit(surv_object ~ NMF_cluster, data = dat.n2)
table(dat.n2$NMF_cluster)
ggsurvplot(fit1, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("SBS7 high", "SBS7 low"),xscale=365.25, break.x.by=730.5,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
?ggsurvplot
#Figure3B
ggsurvplot(fit1, data = dat.n2, pval = FALSE, xlim = c(0,3650), risk.table=TRUE,legend.labs=c("SBS7 high", "SBS7 low"),xscale=365.25, break.x.by=730.5,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_classic() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
fit.coxph <- coxph(surv_object ~  NMF_cluster, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  age_binary, data = dat.n)



#Survival analysis: TCGA dataset
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/03_CLINIC")
dat.n = read.table("survival_excluded.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
colnames(dat.n)

surv_object <- Surv(time = dat.n$futime, event = dat.n$fustat)
fit1 <- survfit(surv_object ~ cluster, data = dat.n)
#Figure3A
ggsurvplot(fit1, data = dat.n, pval = FALSE, xlim = c(0,3650), legend.labs=c("SBS7 high", "SBS7 low"),xscale=365.25, break.x.by=730.5,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
?ggsurvplot
ggsurvplot(fit1, data = dat.n, pval = FALSE, xlim = c(0,3650), risk.table=TRUE,legend.labs=c("SBS7 high", "SBS7 low"),xscale=365.25, break.x.by=730.5,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_classic() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))


###Univariate
fit.coxph <- coxph(surv_object ~  cluster, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  age_binary, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  Mutation_cluster2, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
fit.coxph <- coxph(surv_object ~  sex, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)
dat.n2= dat.n[dat.n$stage!=c("N/A"),]
surv_object2 <- Surv(time = dat.n2$futime, event = dat.n2$fustat)
fit.coxph2 <- coxph(surv_object2 ~  stage, data = dat.n2)
ggforest(fit.coxph2, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)


###Multivariate
###Multivariate
###Multivariate
covariate_names <- c(cluster="cluster",age="age_binary",Mutation_type="Mutation_cluster2",Stage="stage")

map(vars(cluster, age_binary, Mutation_cluster2, stage), function(by)
{
  analyse_multivariate(dat.n2,
                       vars(futime, fustat),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(futime="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor"),
              ggtheme = ggplot2::theme_bw(base_size = 12))

covariate_names <- c(cluster="cluster",age="age_binary",Mutation_type="Mutation_cluster2",Stage="stage",sex="sex")
result=analyse_multivariate(dat.n, vars(futime, fustat),covariates = vars(cluster, age_binary, Mutation_cluster2, stage), covariate_name_dict = covariate_names)
forest_plot(result, factor_labeller = covariate_names,  endpoint_labeller = c(futime="OS"),  orderer = ~order(HR),  labels_displayed = c("endpoint", "factor"),    ggtheme = ggplot2::theme_bw(base_size = 13),      relative_widths = c(1, 1.5, 1))
forest_plot(result, factor_labeller = covariate_names,  endpoint_labeller = c(futime="OS"),  orderer = ~order(HR),  labels_displayed = c("endpoint", "factor"),    ggtheme = ggplot2::theme_bw(base_size = 12))

fit.coxph <- coxph(surv_object ~  age_binary, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)

fit.coxph <- coxph(surv_object ~  Mutation_cluster2, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)

fit.coxph <- coxph(surv_object ~  stage, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)

dat.n2= dat.n[dat.n$stage!=c("N/A"),]
fit.coxph <- coxph(surv_object ~  Mutation_cluster2, data = dat.n)
ggforest(fit.coxph, data = dat.n,fontsize	=1, refLabel="REF",noDigits	=3)


#Survival analysis: TCGA dataset subgroups
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/03_CLINIC")
dat.n = read.table("survival_excluded.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
dat.n2= dat.n[dat.n$stage!=c("N/A"),]
dat.n3= dat.n2[dat.n2$stage==c("Stage I/II"),]
dat.n3= dat.n2[dat.n2$stage==c("Stage III/IV"),]
surv_object6 <- Surv(time = dat.n3$futime, event = dat.n3$fustat)
fit6 <- survfit(surv_object6 ~ cluster, data = dat.n3)
ggsurvplot(fit6, data = dat.n3, pval = FALSE, xlim = c(0,3650), legend.labs=c("SBS7 high", "SBS7 low"),xscale=365.25, break.x.by=730.5,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
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
ggsurvplot(fit6, data = dat.n2, pval = FALSE, xlim = c(0,3650), legend.labs=c("SBS7 high", "SBS7 low"),xscale=365.25, break.x.by=730.5,xlab="Years",palette =c("#0000FF","#FF0000"),conf.int = TRUE, ggtheme =theme_bw() +theme(axis.text.x = element_text(size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_text(size = 14,face="bold"),axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14)))
fit.coxph <- coxph(surv_object6 ~  cluster, data = dat.n2)
ggforest(fit.coxph, data = dat.n2,fontsize	=1, refLabel="REF",noDigits	=3)


##Barplot_proportion
##350/550
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
chisq.test(table(dat.n$Stage_binary, dat.n$cluster,exclude="N/A"))
chisq.test(table(dat.n$Ulcer, dat.n$cluster,exclude="N/A"))
chisq.test(table(dat.n$clark_binary, dat.n$cluster,exclude="N/A"))
dat.n2=dat.n[dat.n$Stage_binary!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$Stage_binary))
#FigureS10
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion") +  theme(legend.title = element_text(size = 14,face="bold"),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")

dat.n2=dat.n[dat.n$Ulcer!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$Ulcer))
dat.n2$Var2=factor(dat.n2$Var2,levels=c("No ulcer","Ulcer"))
#FigureS10
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title = element_text(size = 14,face="bold"),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")

dat.n2=dat.n[dat.n$clark_binary!="N/A",]
dat.n2 <- as.data.frame(table(dat.n2$cluster,dat.n2$clark_binary))
#FigureS10
ggplot(dat.n2, aes(fill=Var2, y=Freq, x=Var1)) +  geom_bar(position="fill", stat="identity",width = 0.7)+theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(color="black",size = 16,face="bold"),axis.text.y = element_text(size = 14), axis.title.y = element_text( size=16, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))  +ylab("Proportion") +  theme(legend.title = element_text(size = 14,face="bold"),legend.text = element_text(size=14), legend.position="top", legend.direction = "horizontal")+labs(fill="")+scale_fill_brewer(palette="Paired")


###Mutation analysis 
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/MAF_Mutect")
MAF_annot=read.table("maf_annotation.csv", sep=",", header=TRUE, row.names =1, stringsAsFactors	=FALSE)
MAF_annot=as.data.frame(MAF_annot)
laml = read.maf(maf = "somatic.maf",clinicalData=MAF_annot)
?read.maf
#Shows sample summry.
getSampleSummary(laml)
summary = getSampleSummary(laml)
##write.csv(x=summary, file="maf_summary.csv", row.names=TRUE)
MAF_annot=read.table("maf_annotation.csv", sep=",", header=TRUE, row.names =1,stringsAsFactors	=FALSE)
?read.table
MAF_annot2=MAF_annot[1:466,]
tail(MAF_annot2)
laml2=subsetMaf(laml, tsb = MAF_annot2$Tumor_Sample_Barcode, query="Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site','Frame_Shift_Del', 'Frame_Shift_Ins')")
getClinicalData(laml2)
laml2@clinical.data
annot_order= as.character(MAF_annot2$order)
fabcolors = c("#FF0000","#0000FF")
names(fabcolors) = c("Non-UV", "UV")
fabcolors = list(FAB_classification = fabcolors)

df2 <- data.frame(supp=rep(c("BRAF", "RAS", "NF1", "Triple WT","KIT"), 2),
                  dose=rep(c("Non-UV", "UV"),each=5),
                  len=c(22.7,13.3,5.3,58.7,12.0,48.2,29.5,9.5,12.8,2.6))
head(df2)
df2$dose=factor(df2$dose,levels=c("UV-low","UV-high"))
df2$supp=factor(df2$supp,levels=c("BRAF", "RAS", "NF1", "Triple WT","KIT"))
#Figure4B
ggplot(data=df2, aes(x=supp, y=len, fill=dose)) + geom_bar(stat="identity", position=position_dodge())+  geom_text(aes(label=len), vjust=1.2, color="white", position = position_dodge(0.9), size=5)+  scale_fill_manual(values=c("#0000FF","#FF0000"))+  theme_classic() +theme(axis.text = element_text( size = 14,face="bold",color="black"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Proportion (%)")+ theme(legend.position = "none") 
##Non-UV 75 BRAF 17 RAS 10 NF1 4 TWT 44 KIT 9
##UV 390 BRAF 188 RAS 115 NF1 37 TWT 50 KIT10

oncoplot(laml2, genes = c("BRAF", "NRAS", "HRAS", "KRAS", "NF1", "CDKN2A", "TP53", "PPP6C", "ARID2", "PTEN", "IDH1", "MAP2K1", "DDX3X", "RAC1", "RB1","KIT"), keepGeneOrder=TRUE,clinicalFeatures=c("order","cluster"),drawRowBar=FALSE, drawColBar=FALSE,removeNonMutated	=FALSE, sortByAnnotation=TRUE,groupAnnotationBySize=FALSE,fill=TRUE,annotationOrder=annot_order, annotationColor=fabcolors,bgCol="#F5F5F5")

#Figure4A
oncoplot(subsetMaf(laml2, tsb = MAF_annot2$Tumor_Sample_Barcode, query="Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site')"), genes = c("BRAF", "NRAS", "HRAS", "KRAS", "NF1", "CDKN2A", "TP53", "PPP6C", "ARID2", "PTEN", "IDH1", "MAP2K1", "DDX3X", "RAC1", "RB1","KIT"), keepGeneOrder=TRUE,clinicalFeatures=c("order","cluster"),drawRowBar=FALSE, drawColBar=FALSE,removeNonMutated	=FALSE, sortByAnnotation=TRUE,groupAnnotationBySize=FALSE,fill=TRUE,annotationOrder=annot_order, annotationColor=fabcolors,bgCol="#F5F5F5")

###copy number alteration analysis

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/")
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
my_comparisons <- list( c("UV-low","UV-high") )
#Figure4C
ggboxplot(dat.n, x = "cluster", y = "NEXUS_Genome",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 16,face="bold",color="black"),axis.text.y = element_text(size = 16), axis.title.y = element_text( size=18, face="bold")) + theme(legend.position = "none") + ylab("Genome with CNA (%)")

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/04_DNA/copynumber")
ak=read.table("SBS7_low_seg.csv", header=TRUE, sep=",")
ak=as.data.frame(ak)
##1000/500
plotFreq(ak,thres.gain=0.2,cex.lab=1.2, cex.axis=1.0)

?plotFreq
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/04_DNA/copynumber")
ak=read.table("SBS7_high_seg.csv", header=TRUE, sep=",")
ak=as.data.frame(ak)
##1000/500
plotFreq(ak,thres.gain=0.2,cex.lab=1.2, cex.axis=1.0)

BiocManager::install("copynumber")
library(copynumber)
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/04_DNA/copynumber")
ak=read.table("total_seg.csv", header=TRUE, sep=",")
ak=as.data.frame(ak)
##1000/500
plotAberration(ak, thres.gain=0.2,sample.labels=FALSE,sep.samples=0)


df2 <- data.frame(supp=rep(c("KIT gain", "TERT gain", "MYC gain", "B2M loss","NOTCH1 loss"), 2),
                  dose=rep(c("Non-UV", "UV"),each=5),
                  len=c(26.7,28.0,58.7,32.0,20.0,11.3,21.2,37.1,15.3,42.5))
head(df2)
df2$dose=factor(df2$dose,levels=c("UV-low","UV-high"))
df2$supp=factor(df2$supp,levels=c("KIT gain", "TERT gain", "MYC gain", "B2M loss","NOTCH1 loss"))
ggplot(data=df2, aes(x=supp, y=len, fill=dose)) + geom_bar(stat="identity", position=position_dodge())+  geom_text(aes(label=len), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5)+  scale_fill_manual(values=c("#FF0000","#0000FF"))+  theme_classic() +theme(axis.text.x = element_text( size = 12,face="bold"),axis.text.y = element_text(size = 12, face="bold"), axis.title.x = element_blank(),axis.title.y = element_text( size=14, face="bold"))+ylab("Proportion of copy alternation (%)")+ theme(legend.position = "none") 


##Non-UV 75 BRAF 17 RAS 10 NF1 4 TWT 44 KIT 9
##UV 390 BRAF 188 RAS 115 NF1 37 TWT 50 KIT10


###DEG analysis 
###DEG analysis 
###DEG analysis 


setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/04_RNA/DEG")
x=read.table("filtered_RNA_raw.txt", sep="\t", row.names=1, header=TRUE)
group <- read.table("rna_annotation.txt", sep="\t", row.names=1, header=FALSE)
group$V2=factor(group$V2,levels=c("Non-UV", "UV"))
y <- DGEList(counts=x, group=group$V2)
countsPerMillion <- cpm(y)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 117.25)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+group, data=y$samples)
design
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
ggplot(input, aes(logFC, -log10(FDR))) + geom_point(aes(col=significant)) +  scale_color_manual(values=c("black", "red")) + theme(text = element_text(size=12))+geom_label_repel(data=head(input,15), aes(label=gene))+ theme_classic() +theme(axis.text.x = element_text(size=14,face="bold"),axis.text.y = element_text(size=14, face="bold"), axis.title = element_text(size=16, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14))+ theme(legend.position = "none")

###Cibersortx heatmap
###Cibersortx heatmap

dat.n=read.table("cibersortx.csv", sep=",", quote="",row.names=1, header=TRUE)
dat.n
ggplot(dat.n, aes(x=factor(1:length(Absolute)), y=Absolute)) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ylab("Absolute score from Cibersortx")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

###Cibersortx test 
###Cibersortx test
###Cibersortx test
dat.n=read.table("cibersortx.csv", sep=",", quote="",row.names=1, header=TRUE)
dat.n2=as.matrix(t(dat.n)[1:22,])
p_values <- vector("list", nrow(dat.n2))
for(i in seq_along(1: nrow(dat.n2))){
  p_values[i] = wilcox.test(dat.n2[i,1:75],dat.n2[i,76:465], paired = FALSE, alternative = "two.sided")$p.value
}
p_values = data.frame(p_values = sapply(p_values, c))
rownames(p_values)=rownames(dat.n2)

dat.n=read.table("cibersortx_relative.csv", sep=",", quote="",row.names=1, header=TRUE)
dat.n2=as.matrix(t(dat.n)[1:22,])
p_values <- vector("list", nrow(dat.n2))
for(i in seq_along(1: nrow(dat.n2))){
  p_values[i] = wilcox.test(dat.n2[i,1:75],dat.n2[i,76:465], paired = FALSE, alternative = "two.sided")$p.value
}

p_values = data.frame(p_values = sapply(p_values, c))
rownames(p_values)=rownames(dat.n2)



##Cibersort visualization 
#Figure S11 
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/04_RNA/cibersort")
dat.n=read.table("cibersortx.csv", sep=",", quote="",row.names=1, header=TRUE)
##dat.n=read.table("cibersortx_relative.csv", sep=",", quote="",row.names=1, header=TRUE)
dat.n$cluster=rep(c("UV-low", "UV-high"), c(75,390))
table(dat.n$cluster)
my_comparisons <- list( c("UV-low", "UV-high"))
dat.n=as.data.frame(dat.n)

dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))
dat.n
ggboxplot(dat.n, x = "cluster", y = "Absolute",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")+ylab("Total immune cell fraction")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

ggboxplot(dat.n, x = "cluster", y = "Macrophages_M1",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

ggboxplot(dat.n, x = "cluster", y = "T_cells_CD4_memory_activated",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

ggboxplot(dat.n, x = "cluster", y = "T_cells_regulatory",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

ggboxplot(dat.n, x = "cluster", y = "T_cells_CD8",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text( size = 16,face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ theme(legend.position="none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

###Tumor immunity markers
#Figure 4F
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("NEW_cluster_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels= c("UV-low","UV-high"))

ggboxplot(dat.n, x = "cluster", y = "purity",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none") +ylab("Tumor purity")+ theme(legend.position = "none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

ggboxplot(dat.n, x = "cluster", y = "Leukocyte",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA) + theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+ylab("Leukocyte fraction")+ theme(legend.position = "none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

dat.n2=dat.n[is.na(dat.n$CYT_score_RSEM)!=TRUE,]
dat.n2$CYT_score_RSEM=as.numeric(as.character(dat.n2$CYT_score_RSEM))
ggboxplot(dat.n2, x = "cluster", y = "CYT_score_RSEM",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA) + theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+ylab("CYT score (RSEM)")+ theme(legend.position = "none")+scale_y_continuous(trans='log10')
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")


dat.n2=dat.n[is.na(dat.n$TIL_percentage)!=TRUE,]
dat.n2$TIL_percentage=as.numeric(as.character(dat.n2$TIL_percentage))
ggboxplot(dat.n2, x = "cluster", y = "TIL_percentage",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+ylab("Proportion of tumor infiltrating leukocytes")+theme(legend.position = "none")
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")

dat.n2=dat.n[is.na(dat.n$TCR)!=TRUE,]
dat.n2$TCR=as.numeric(as.character(dat.n2$TCR))
ggboxplot(dat.n2, x = "cluster", y = "TCR",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ theme_classic() +theme(axis.text.x = element_text(size = 16, face="bold",color="black"),axis.text.y = element_text(size = 14), axis.title.x = element_blank(),axis.title.y = element_text(size=16, face="bold"))+  theme(legend.position = "none")+  theme(legend.position = "none") +ylab("TCR richness")+ theme(legend.position = "none")+scale_y_continuous(trans='log10')
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")



dat.n2=dat.n[dat.n$PDL1!="N/A",]
dat.n2$PDL1=as.numeric(as.character(dat.n2$PDL1))

ggboxplot(dat.n2, x = "cluster", y = "PDL1",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test") + theme_classic() +theme(axis.text.x = element_text(size = 12, face="bold"),axis.text.y = element_text(size = 12), axis.title.x = element_blank(),axis.title.y = element_text(size=12, face="bold"))+ylab("PD-L1 normalized expression (RPPA)")+ theme(legend.position = "none")
?ggbarstats


###panel sequencing data preparation
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
mutation=read.table("MSK-SKCM_filter_SNP.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
head(sample)

list=sample$file
for (list in list) {
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/SNP")
  vcf.obj <- FIREVAT::ReadVCF(list,genome="hg19")
  target.mut.sigs <- c("SBS7a")
  mutalisk.results <- RunMutalisk(vcf.obj,df.ref.mut.sigs = GetPCAWGMutSigs(),target.mut.sigs =target.mut.sigs ,method="all",verbose=TRUE)
  ?RunMutalisk
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/TRI")
  write.table(x=t(mutalisk.results$sub.types.spectrum), file=list, sep=" ")
}
sig=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
sig=as.data.frame(sig)
colnames(sig)=colnames(sigs.input)
sample=rownames(sig)

for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = list1)
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/COSMIC_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/COSMIC_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, sample.id = "Sample", chr = "chr",pos = "pos",ref = "ref",alt = "alt")
sig=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
head(sig)
sig=as.data.frame(sig)
colnames(sig)=colnames(sigs.input)
sample=rownames(sig)

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
sig2 = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(sigs.input)

for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sig, signatures.ref =sig2, contexts.needed = FALSE, associated	=c("Basis1", "Basis2"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/basis_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/basis_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}


setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
###Figure 5A
sig2 = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(sigs.input)
rownames(sig2)
dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = "MSK_SBS7_low")


setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
pdf('MSK_SBS7_low.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()

###Figure 5A
dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id = "MSK_SBS7_high")

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
pdf('MSK_SBS7_high.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()



###Panel sequencing (MSK cohort) Visualization 
###Figure 5C
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
colnames(dat.n)
dat.n2=t(as.data.frame(dat.n[,6:8]))
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=rownames(dat.n)
colnames(my_sample_col)="Cluster"
dat.n2=dat.n2[,order(sample$SigB)]
rownames(dat.n2)=c("SigA","SigB","Unknown")
my_gene_col <- data.frame(sample = rep(c("SigA","SigB","Unknown"), c(1,1,1)))
rownames(my_gene_col)=rownames(dat.n2)
colnames(my_gene_col)="Basis"
my_gene_col$Basis=factor(my_gene_col$Basis, levels= c("SigA","SigB","Unknown"))
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Cluster=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_names_row=FALSE,color = colorRampPalette(brewer.pal(n = 11, name ="YlOrRd"))(100),annotation_col =my_sample_col,annotation_names_col=FALSE,annotation_colors=ann_colors,legend=F,annotation_legend = F)



###Panel sequencing (MSK cohort) Visualization 
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(dat.n)
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Cluster"
dat.n=dat.n[,order(sample$Basis2)]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Type = Var1, Cluster=Var2)
##1200/600
#Figure S13
pheatmap(dat.n,show_rownames=FALSE,show_colnames=FALSE,cluster_cols=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_row =my_gene_col,annotation_col =my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)

dat.n=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
dat.n2=dat.n[order(sample$Basis2),]
sum=cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% tidyr::gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")


###Basis1,2 signature
###Basis1,2 signature
###Basis1,2 signature
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n=as.data.frame(dat.n[,6:8])
dat.n=t(dat.n)
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Cluster"


dat.n=dat.n[,order(sample$Basis2)]
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list(Cluster=Var2)

##1200/600
#Figure 5C, upper 
pheatmap(dat.n,show_rownames=TRUE,show_colnames=FALSE,cluster_cols=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col =my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)

###COSMIC signature
###COSMIC signature
###COSMIC signature
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n=t(as.data.frame(dat.n[,12:31]))
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
my_sample_col=as.data.frame(sample$cluster)
rownames(my_sample_col)=colnames(dat.n)
colnames(my_sample_col)="Cluster"
dat.n=dat.n[,order(sample$SigB)]
Var2 = c("#0000FF","#FF0000")
names(Var2) = c("Non-UV", "UV")
ann_colors = list( Cluster=Var2)
##1200/600
#Figure 5C, lower
pheatmap(dat.n,show_rownames=TRUE,show_colnames=FALSE,cluster_cols=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col =my_sample_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,fontsize = 12,legend=F,annotation_legend = F)


dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(dat.n$Basis2),]
colnames(dat.n2)
ggplot(dat.n2, aes(x=factor(1:length(mutation)), y=as.numeric(as.character(mutation)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ scale_y_continuous(trans='log10',sec.axis = dup_axis()) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

ggplot(dat.n2, aes(x=factor(1:length(basis_cos)), y=as.numeric(as.character(basis_cos)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

ggplot(dat.n2, aes(x=factor(1:length(COSMIC_cos)), y=as.numeric(as.character(COSMIC_cos)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 


##400/600
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("panel_sig_annotation.txt", sep="\t", row.names=1, header=TRUE)
dat.n$cluster=factor(dat.n$cluster, levels=c("UV-low","UV-high"))
colnames(dat.n)
#Figure5D
ggboxplot(dat.n, x = "cluster", y = "mutation",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text = element_text( size = 14,face="bold",color="black"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Mutation counts")+ theme(legend.position = "none") +scale_y_continuous(trans='log10')
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")
#Figure5D
ggboxplot(dat.n, x = "cluster", y = "Sum_SBS7",color = "cluster", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "cluster",outlier.shape = NA)+  theme_classic() +theme(axis.text = element_text( size = 14,face="bold",color="black"), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"))+ylab("Sum of SBS7")+ theme(legend.position = "none") 
+ stat_compare_means(comparisons = my_comparisons, aes(label = ..p.format..), method ="wilcox.test")



dat.n=read.table("panel_sig_tri.txt", sep="\t", row.names=1, header=TRUE)
sample=read.table("panel_sig_annotation.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
dat.n2=dat.n[order(sample$Basis2),]
dim(dat.n2)
head(dat.n2)
sum=rbind(colSums(dat.n2[1:12,])/12,colSums(dat.n2[13:245,])/233)
write.table(x=sum, file="average.txt",quote = FALSE, sep="\t",row.names = TRUE,col.names = TRUE)

###PANEL simulation data preparation
dat.n=read.table("simulate_deconstruc.txt", sep="\t",  header=TRUE)
sim.sig <- mut.to.sigs.input(mut.ref = dat.n, sample.id = "Sample", chr = "chr",pos = "pos", 
                             ref = "ref", alt = "alt",bsg = BSgenome.Hsapiens.UCSC.hg38)

sample=dat.n$Sample

for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sim.sig, signatures.ref =signatures.exome.cosmic.v3.may2019, contexts.needed = TRUE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"),  sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/SIM_COSMIC_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/SIM_COSMIC_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/SIM_COSMIC_TRI")
  write.table(x=tumor, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}


setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
sig2 = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(sigs.input)


for (list1 in sample) {
  dec = whichSignatures(tumor.ref = sim.sig, signatures.ref =sig2, contexts.needed = TRUE, associated	=c("Basis1", "Basis2"), sample.id = list1)
  dec$weights
  tumor=dec$tumor
  deconstruct=dec$product
  cos=cos.sim(tumor,deconstruct)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/SIM_basis_cos")
  write.table(x=cos, file=list1, quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
  setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL/SIM_basis_sig")
  weight=dec$weights
  write.table(x=weight, file=list1,quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE)}

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("simulate_tri2.txt", sep="\t", row.names=1, header=TRUE)
sample=read.table("simulate_annot2.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
head(sample)
dat.n2=dat.n[order(sample$orig.order),]
table(sample$cluster)
table(sample$orig.ident)
sum=rbind(colSums(dat.n2[1:54,])/54,colSums(dat.n2[55:383,])/329)
write.table(x=sum, file="average.txt",quote = FALSE, sep="\t",row.names = TRUE,col.names = TRUE)

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("simulate_tri2.txt", sep="\t", row.names=1, header=TRUE)
sample=read.table("simulate_annot2.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
head(sample)
dat.n2=dat.n[order(sample$Basis2),]
head()
sum=rbind(colSums(dat.n2[1:13,])/13,colSums(dat.n2[14:383,])/370)
write.table(x=sum, file="average.txt",quote = FALSE, sep="\t",row.names = TRUE,col.names = TRUE)



#Comparison between WES and panel simulation of WES
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("simulate_annot2.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
dat.n=dat.n[order(dat.n$orig.order),]
head(dat.n)
tail(dat.n)
colnames(dat.n)
median(dat.n$mut.retention)
#Figure S12
ggplot(dat.n, aes(x=factor(1:length(mut.retention)), y=as.numeric(as.character(mut.retention)))) +geom_bar(stat = "identity", colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Mutation counts")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))

setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
dat.n=read.table("simulate_tri2.txt", sep="\t", row.names=1, header=TRUE)
dat.n2=dat.n[order(sample$orig.order),]
setwd("/home/kysbbubbu/skcm_sig_MAF")
dat.n=read.table("MAF_heatmap.csv", sep=",", row.names=1, header=TRUE)
dat.n3=dat.n[rownames(dat.n2),]
dat.n2=as.matrix(dat.n2)
dat.n3=as.matrix(dat.n3)
cos=cos.sim(dat.n2,dat.n3)
cos.matrix=diag(cos)
cos.matrix=as.data.frame(cos.matrix)
min(cos.matrix)
colnames(cos.matrix)="cos.sim"
median(cos.matrix$cos.sim)
#Figure S12
ggplot(cos.matrix, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="solid",colour="black"))+scale_y_continuous(limits=c(0,1),breaks = scales::pretty_breaks(n = 2)) 

dat.n2=as.data.frame(t(dat.n2))
dat.n3=as.data.frame(t(dat.n3))
my_sample_col=as.data.frame(colnames(dat.n2))
colnames(my_sample_col)="Cluster"
rownames(my_sample_col)=colnames(dat.n2)
my_sample_col$Cluster=rep(c("SBS7 low","SBS7 high"), c(13,370))
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n2)
colnames(my_gene_col)="Type"
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("SBS7 low", "SBS7 high")
ann_colors = list(Type = Var1, Cluster=Var2)
#Figure S12
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend=TRUE)
#Figure S12
pheatmap(dat.n3,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,legend=FALSE,annotation_legend=TRUE)

dat.n2=t(dat.n2)
sum=cbind(rowSums(dat.n2[,1:16]),rowSums(dat.n2[,17:32]),rowSums(dat.n2[,33:48]),rowSums(dat.n2[,49:64]),rowSums(dat.n2[,65:80]),rowSums(dat.n2[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% tidyr::gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
#Figure S12
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")
dat.n3=t(dat.n3)
sum=cbind(rowSums(dat.n3[,1:16]),rowSums(dat.n3[,17:32]),rowSums(dat.n3[,33:48]),rowSums(dat.n3[,49:64]),rowSums(dat.n3[,65:80]),rowSums(dat.n3[,81:96]))
colnames(sum)=c("C.A","C.G","C.T","T.A","T.C","T.G")
sum=as.data.frame(sum)
datm <- sum %>% dplyr::mutate(ind = factor(row_number())) %>% tidyr::gather(variable, value, -ind)
datm$variable  <- plyr::mapvalues(x = datm$variable, from = c("C.A","C.G","C.T","T.A","T.C","T.G"), to = c("C>A","C>G","C>T","T>A","T>C","T>G"))
##650/150
#Figure S12
ggplot(datm, aes(x=ind, y=value, fill=variable),space=0) +geom_bar(stat = "identity", position = "fill",width=1) + theme(text = element_blank()) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),plot.background=element_blank(),panel.border=element_blank())+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + theme(legend.position = "none")


#Visualization  of  panel sequencing simulation data
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/05_VALIDATION")
sig2 = read.table("deconstruc.csv", row.names=1,sep=",",header = T, stringsAsFactors =FALSE)
sig2=as.data.frame(t(sig2))
colnames(sig2)=colnames(sigs.input)
rownames(sig2)
dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id =  "NEW_SIM_SBS7_low" )
#Figure 5B
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
pdf('NEW_SIM_SBS7_low.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()
dec = whichSignatures(tumor.ref = sig2, signatures.ref =signatures.exome.cosmic.v3.may2019,  contexts.needed = FALSE, associated	=c("SBS1", "SBS2","SBS3", "SBS5", "SBS7a","SBS7b","SBS7c","SBS7d","SBS9", "SBS11","SBS13","SBS14", "SBS17a","SBS17b","SBS31","SBS36", "SBS38","SBS40","SBS43", "SBS45","SBS49","SBS52","SBS55","SBS58"), sample.id =  "NEW_SIM_SBS7_high" )
#Figure 5B
setwd("/data/MRC1_data4/kysbbubbu/skcm_sig_MAF/06_PANEL")
pdf('NEW_SIM_SBS7_high.pdf', width = 10, height = 10)
chart<-plotSignatures(dec)
dev.off()
