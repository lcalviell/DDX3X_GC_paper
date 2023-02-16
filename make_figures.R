.libPaths(c(.libPaths(),"/home/lorenzo.calviello/R/x86_64-pc-linux-gnu-library/4.0/"))

library(RiboseQC)
library(INSPEcT)
library(cowplot)
library(corrplot)
library(corrgram)


setwd("/group/calviello/Lorenzo/DDX3_decay/To_Github/")




source("data/new_riboseqc_all_HT.R")
annotation_file<-"data/gencode.v32.annotation.gtf.gz_Rannot"
load(annotation_file)

inspres<-readRDS("data/INSPEcT_res")
gns<-GTF_annotation$trann$gene_name[match(rownames(inspres@ratePvals),GTF_annotation$trann$gene_id)]
#plotGene(inspres, ix=featureNames(inspres)[gns=="TLE4"],priors = F,constantModel = T)
summary_multiDE<-get(load("data/degron_multiDE_results_summary.RData"))

gcqnt_orig<-summary_multiDE$uniq[[1]][,c("gene_id","gene_name","GCpct_cds")]


trans <- function(x) -log(x, 10)
inv <- function(x) 10^(-x)

df<-summary_multiDE$uniq


deg1<-viewModelRates(inspres,"synthesis")
deg1cp<-deg1
deg1<-t(apply(deg1,1,function(x){log2(x/x[1])}))
deg2<-viewModelRates(inspres,"processing")
deg2cp<-deg2
deg2<-t(apply(deg2,1,function(x){log2(x/x[1])}))
deg3<-viewModelRates(inspres,"degradation")
deg3cp<-deg3
deg3<-t(apply(deg3,1,function(x){log2(x/x[1])}))

deg3[abs(deg3)>8]<-NA
deg2[abs(deg2)>8]<-NA
deg1[abs(deg1)>8]<-NA

deg3<-deg3[rownames(deg3)%in%df[[1]]$gene_id,]
deg2<-deg2[rownames(deg2)%in%df[[1]]$gene_id,]
deg1<-deg1[rownames(deg1)%in%df[[1]]$gene_id,]

deg3cp<-deg3cp[rownames(deg3cp)%in%df[[1]]$gene_id,]
deg2cp<-deg2cp[rownames(deg2cp)%in%df[[1]]$gene_id,]
deg1cp<-deg1cp[rownames(deg1cp)%in%df[[1]]$gene_id,]

for(zz in 1:length(df)){
  dff<-df[[zz]]
  #mtcc<-match(dff$gene_id,rownames(deg3))
  mtcc<-match(rownames(deg3),dff$gene_id)
  okk<-which(!is.na(mtcc))
  dff$base_deg<-NA
  dff$base_syn<-NA
  dff$base_proc<-NA
  dff$base_deg[mtcc]<-deg3cp[,1]
  dff$log2_deg<-NA
  dff$base_syn[mtcc]<-deg1cp[,1]
  dff$log2_syn<-NA
  dff$base_proc[mtcc]<-deg2cp[,1]
  dff$log2_proc<-NA
  
  dff$log2_deg[mtcc]<-deg3[,zz+1]
  dff$log2_syn[mtcc]<-deg1[,zz+1]
  dff$log2_proc[mtcc]<-deg2[,zz+1]
  
  
  df[[zz]]<-dff
}

sig_genes<-apply(inspres@ratePvals,1,function(x){min(x,na.rm = T)})
sig_genes<-names(sig_genes[which(sig_genes<.05)])
trans <- function(x) {reso=-log(x, 10);reso[is.na(reso)]=0;return(reso)}
inv <- function(x) 10^(-x)

df<-do.call(df,what = rbind)
rownames(df)<-NULL
df$baseline_RiboRNA[is.infinite(df$baseline_RiboRNA)]<-NA
i="baseline_RiboRNA"
df[,i]<-log2(df[,i])
df$baseline_RiboRNA[is.infinite(df$baseline_RiboRNA)]<-NA

df[,i]<-(df[,i]-mean(df[complete.cases(df[,i]),i]))/sd(df[complete.cases(df[,i]),i])

df$single_exon<-df$n_exons==1
df$IntronExon_log2FC[df$intronlen==0]<-NA
df$IntronExon_padj[df$intronlen==0]<-NA

trans <- function(x) -log(x, 10)
inv <- function(x) 10^(-x)
dfa<-df
if(is.null(dfa$TPM_RNA)){
  baseMean_len<-dfa$RNA_baseMean/(dfa$exlen)
  baseMean_len<-baseMean_len*(1000000/sum(baseMean_len,na.rm = T))
  dfa$TPM_RNA<-baseMean_len
}

if(is.null(dfa$TPM_Ribo)){
  baseMean_len<-dfa$Ribo_baseMean/(dfa$exlen)
  baseMean_len<-baseMean_len*(1000000/sum(baseMean_len,na.rm = T))
  dfa$TPM_Ribo<-baseMean_len
}
dfa$baseline_RiboRNA<-log(dfa$TPM_Ribo/dfa$TPM_RNA)
dfa$baseline_RiboRNA[is.infinite(dfa$baseline_RiboRNA)]<-NA

dfa[,"baseline_RiboRNA"]<-(dfa[,"baseline_RiboRNA"]-mean(dfa[complete.cases(dfa[,"baseline_RiboRNA"]),"baseline_RiboRNA"]))/sd(dfa[complete.cases(dfa[,"baseline_RiboRNA"]),"baseline_RiboRNA"])

dfa$baseline_IntronExon<-log(dfa$TPM_Intron/dfa$TPM_RNA)
dfa$baseline_IntronExon[is.infinite(dfa$baseline_IntronExon)]<-NA

dfa[,"baseline_IntronExon"]<-(dfa[,"baseline_IntronExon"]-mean(dfa[complete.cases(dfa[,"baseline_IntronExon"]),"baseline_IntronExon"]))/sd(dfa[complete.cases(dfa[,"baseline_IntronExon"]),"baseline_IntronExon"])


dfa<-dfa[which(dfa$TPM_RNA>3 & dfa$TPM_Ribo>3),]

colsi<-alpha(c("blue","dark red","gray43","gray10","azure2","dark gray"),c(.8,.8,.8,.8,.3,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed","unchanging")

dfa$tx_type_RiboRNA<-"mixed"
dfa$tx_type_IntronExon<-dfa$tx_type_RiboRNA
dfa$tx_type_RiboRNA[which(dfa$Ribo_padj>.1 & dfa$RNA_padj>.1 | is.na(dfa$RNA_padj) | is.na(dfa$Ribo_padj) )]<-"unchanging"
dfa$tx_type_RiboRNA[dfa$RNA_padj<.01 & dfa$RNA_log2FC<0 & dfa$Ribo_log2FC<0]<-"Concordant_down"
dfa$tx_type_RiboRNA[dfa$RNA_padj<.01 & dfa$RNA_log2FC>0 & dfa$Ribo_log2FC>0]<-"Concordant_up"

dfa$tx_type_IntronExon[dfa$RNA_padj<.01 & dfa$RNA_log2FC<0 & dfa$Intron_log2FC<0]<-"Concordant_down"
dfa$tx_type_IntronExon[dfa$RNA_padj<.01 & dfa$RNA_log2FC>0 & dfa$Intron_log2FC>0]<-"Concordant_up"
dfa$tx_type_IntronExon[which(dfa$RNA_padj>.1 & dfa$Intron_padj>.1 | is.na(dfa$RNA_padj) | is.na(dfa$Intron_padj) )]<-"unchanging"

dfa$tx_type_RiboRNA[dfa$RiboRNA_padj<.05 & dfa$RiboRNA_log2FC<0]<-"TE_down"
dfa$tx_type_RiboRNA[dfa$RiboRNA_padj<.05 & dfa$RiboRNA_log2FC>0]<-"TE_up"
dfa$tx_type_RiboRNA<-factor(dfa$tx_type_RiboRNA,levels=c("TE_down","TE_up","Concordant_down","Concordant_up","mixed","unchanging"))

dfa$tx_type_IntronExon[which(dfa$IntronExon_padj>.1 & dfa$RNA_padj>.1)]<-"unchanging"
dfa$tx_type_IntronExon[dfa$IntronExon_padj<.05 & dfa$IntronExon_log2FC<0]<-"premRNA_down"
dfa$tx_type_IntronExon[dfa$IntronExon_padj<.05 & dfa$IntronExon_log2FC>0]<-"premRNA_up"
dfa$tx_type_IntronExon<-factor(dfa$tx_type_IntronExon,levels=c("premRNA_down","premRNA_up","Concordant_down","Concordant_up","mixed_NS"))
livv<-levels(dfa$tx_type_RiboRNA)
livv<-livv[!livv%in%c("mixed","unchanging")]
okgns<-c()
for(li in livv){
  dfo<-dfa[dfa$tx_type_RiboRNA==li,]
  dfo$RNA_padj[is.na(dfo$RNA_padj)]<-1
  dfo$RiboRNA_padj[is.na(dfo$RiboRNA_padj)]<-1
  if(grepl(li,pattern = "Concor")){okgns<-c(okgns,dfo[order(dfo$RNA_padj,decreasing = F),"gene_name"][1:2])}
  if(grepl(li,pattern = "TE_")){okgns<-c(okgns,dfo[order(dfo$RiboRNA_padj,decreasing = F),"gene_name"][1:2])}
}
dfaaa<-split(dfa,dfa$tx_type_RiboRNA)
dfaaa[["mixed"]]<-NULL
dfaaa[["unchanging"]]<-NULL

dfa$gene_name_ok<-dfa$gene_name
dfa$gene_name<-NA
dfa$gene_name[dfa$gene_name_ok%in%okgns]<-dfa$gene_name_ok[dfa$gene_name_ok%in%okgns]

dfa$RiboRNA_padj[is.na(dfa$RiboRNA_padj)]<-1
dfa$RiboRNA_padj[dfa$RiboRNA_padj<10e-50]<-10e-50

dfa$IntronExon_padj[is.na(dfa$IntronExon_padj)]<-1
dfa$IntronExon_padj[dfa$IntronExon_padj<10e-50]<-10e-50

dfa$RNA_padj[is.na(dfa$RNA_padj)]<-1
dfa$RNA_padj[dfa$RNA_padj<10e-50]<-10e-50
dfa$experiment<-gsub(dfa$experiment,pattern = "_vs_DMSO",replacement = "")
axm<-max(abs(c(dfa$RNA_log2FC,dfa$Ribo_log2FC)),na.rm = T)

colsi<-alpha(c("blue","dark red","gray43","gray10","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_NS")
dfaz<-dfa
dfaz$tx_type_RiboRNA<-as.character(dfaz$tx_type_RiboRNA)
dfaz$tx_type_RiboRNA[dfaz$tx_type_RiboRNA%in%c("mixed","unchanging")]<-"mixed_NS"
dfaz$tx_type_RiboRNA<-factor(dfaz$tx_type_RiboRNA,levels=names(colsi))
dfaz$tx_type_RiboRNA<-factor(dfaz$tx_type_RiboRNA,levels=names(colsi))
dfaz$tx_type_RiboRNA[dfaz$tx_type_RiboRNA%in%c("Concordant_down","Concordant_up")]<-"mixed_NS"
a<-ggplot(dfaz,aes(x=RNA_log2FC,y=Ribo_log2FC,color=tx_type_RiboRNA,size=RiboRNA_padj,label=gene_name,gene_symbol=gene_name_ok)) + geom_point()
if(sum(dfa$gene_biotype%in%c("lncRNA","ncRNA"))>500){
  a<-ggplot(dfa,aes(x=RNA_log2FC,y=Ribo_log2FC,color=tx_type_RiboRNA,size=RiboRNA_padj,label=gene_name,gene_symbol=gene_name_ok)) + geom_point()
}
a<-a + theme_bw() +
  ylab(paste("Ribo-seq "," log2FC",sep = "")) +
  xlab(paste("RNA-seq log2FC",sep = "")) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"Tx_class") +
  scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv,
                                                                               domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))
axm=5
a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + coord_cartesian(xlim = c(-axm,axm),ylim=c(-axm,axm))
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1)
pdf(file = "Figures/Fig1A.pdf",width = 22,height = 5)
ariborna
dev.off()

dfaz$time<-as.integer(gsub(gsub(dfaz$experiment,pattern = "time",replacement = ""),pattern = "h",replacement = ""))

dfaz2<-dfaz[dfaz$tx_type_RiboRNA%in%c("TE_up","TE_down"),]
diffaz<-aggregate(dfaz2$RiboRNA_log2FC,by=list(dfaz2$tx_type_RiboRNA,dfaz2$time),mean)
colnames(diffaz)<-c("tx_type_RiboRNA","time","RiboRNA_log2FC")
diffaz$n_sig<-c(t(table(dfaz2$time,dfaz2$tx_type_RiboRNA)[,1:2]))
a<-ggplot(diffaz,aes(x=time,y=RiboRNA_log2FC,color=tx_type_RiboRNA,size=n_sig)) + geom_point()
a<-a + theme_classic() +
  ylab("delta TE\n(average)") +
  xlab("time after degron induction") +
  geom_hline(yintercept =c(-1,0,1),linetype=2 ) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"Tx_class") +
  scale_size_continuous(range = c(2,15),breaks =c(50,200,500),"number of\nregulated genes")

ariborna<-a + coord_cartesian(ylim=c(-1.7,1.7))
pdf(file = "Figures/Fig1B.pdf",width = 10,height = 5)
ariborna  +  scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48))) + geom_line(size=1)
dev.off()

dfa2<-dfa
gnsok<-unique(dfa2$gene_id[which(dfa$RNA_padj<.01)])
#dfa2<-dfa2[dfa2$gene_id%in%gnsok,]
dfa2$time<-as.integer(gsub(gsub(dfa2$experiment,pattern = "time",replacement = ""),pattern = "h",replacement = ""))
df<-data.frame(fcs=c(dfa2$log2_syn,dfa2$log2_proc,dfa2$log2_deg,dfa2$Ribo_log2FC,dfa2$RNA_log2FC))
df$gene_id<-rep(dfa2$gene_id,5)
dfa2$GCpct_cds<-dfa2$GCpct_cds*100

livo<-setNames(nm = dfa2[dfa2$experiment=="48h",]$gene_id,dfa2[dfa2$experiment=="48h",]$tx_type_RiboRNA)

#livo<-setNames(nm = dfa[dfa$experiment=="time48h",]$gene_id,interaction(dfa[dfa$experiment=="time48h",]$xy_RiboRNA,dfa[dfa$experiment=="time48h",]$Ribo_log2FC>0))
#livo<-setNames(nm = dfa2[dfa2$experiment=="time48h",]$gene_id,dfa2[dfa2$experiment=="time48h",]$RiboRNA_log2FC>0)

df$tx_type_RiboRNA<-livo[df$gene_id]
df$time<-rep(dfa2$time,5)
df$assay<-rep(c("synthesis","processing","degradation","translation","abundance"),each=length(dfa2[,1]))

gcqnt<-setNames(cut(dfa2$GCpct_cds[dfa2$time==4],breaks = quantile(dfa2$GCpct_cds[dfa2$time==4],prob=c(seq(0,1,length.out = 6))),include.lowest = T,right = T),nm = dfa2$gene_id[dfa2$time==4])
gcqnt<-gsub(gsub(gsub(gsub(as.character(gcqnt),pattern = "[(]",replacement = ""),pattern = "\\[",replacement = ""),pattern = "\\]",replacement = ""),pattern = ",",replacement = "-")
gcqnt<-factor(gcqnt)
gcqnt<-setNames(gcqnt,nm = dfa2$gene_id[dfa2$time==4])
df$gcqnt<-gcqnt[df$gene_id]

#adjust here to have all combinations for plotting at zero
aa<-unique(df[which(df$fcs==0),])[rep(1:5,each= 30),]
aa$assay<-rep(c("synthesis","processing","degradation","translation","abundance"),each=30)
aa$time<-0
aa$tx_type_RiboRNA<-rep(levels(aa$tx_type_RiboRNA),25)
aa$gcqnt<-rep(levels(aa$gcqnt),30)
df<-rbind(df,aa)

df$assay<-factor(df$assay,levels=c("abundance","translation","synthesis","processing","degradation"))


df<-df[df$gene_id%in%sig_genes | df$time==0,]
dfz<-df
dfz<-dfz[dfz$tx_type_RiboRNA%in%c("TE_up","TE_down","unchanging"),]
colsi<-alpha(c("blue","dark red","dark gray"),c(.8,.8,.8))
names(colsi)<-c("TE_down","TE_up","unchanging")

dfz<-dfz[dfz$assay%in%c("abundance","translation"),]
dfz$assay<-as.character(dfz$assay)

dfz$assay[dfz$assay=="abundance"]<-"RNA-seq"
dfz$assay[dfz$assay=="translation"]<-"Ribo-seq"

dfz$assay<-factor(dfz$assay,levels=c("Ribo-seq","RNA-seq"))
dfz$inter<-interaction(dfz$tx_type_RiboRNA,dfz$assay)
a<-ggplot(dfz,aes(x=time,y=fcs,color=tx_type_RiboRNA,group=inter,linetype=assay,shape=assay))  + stat_summary(geom="line",size=2) +stat_summary(size=2)
a<-a + theme_classic() +
  ylab("log2FC") +
  xlab("hours post degron induction") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"Tx_class") + scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48)))
#a<-a + facet_wrap(.~assay)
a<-a + geom_hline(yintercept = 0,size=.5,lty=2)

pdf(file = "Figures/Fig1C.pdf",width = 10,height = 5)
a
dev.off()

colsi<-alpha(colorRampPalette(colors = c("forestgreen","grey","purple"))(5),c(.8,.8,.8,.8,.8))
#colsi<-alpha(c("blue","dark red","dark gray"),c(.8,.8,.8,.8,.3))
#names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_NS")

df<-df[df$assay!="translation",]
df<-df[df$assay!="abundance",]

df$assay<-factor(df$assay,levels=c("synthesis","processing","degradation"))


df<-df[df$tx_type_RiboRNA%in%c("TE_up","TE_down","unchanging"),]
colsi<-alpha(c("blue","dark red","dark gray"),c(.8,.8,.8))
names(colsi)<-c("TE_down","TE_up","unchanging")

a<-ggplot(df,aes(x=time,y=fcs,color=gcqnt,group=gcqnt))  + stat_summary(geom="line") +stat_summary()

colsi<-alpha(colorRampPalette(colors = c("forestgreen","grey","purple"))(5),c(.8,.8,.8,.8,.8))


a<-a + theme_classic() +
  ylab("log2FC") +
  xlab("hours post degron induction") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48))) +
  scale_color_manual(values = colsi,"%GC in CDS")

a<-a + facet_wrap(.~assay,nrow = 1)
a<-a + geom_hline(yintercept = 0,size=.5,lty=2)
pdf(file = "Figures/Supplementary_Figure_5.pdf",width = 12,height = 3.5)
a
dev.off()

dfabc<-df
dfabc<-dfabc[dfabc$assay=="degradation",]
a<-ggplot(dfabc,aes(x=time,y=fcs,color=gcqnt,group=gcqnt))  + stat_summary(geom="line") +stat_summary()

colsi<-alpha(colorRampPalette(colors = c("forestgreen","grey","purple"))(5),c(.8,.8,.8,.8,.8))


a<-a + theme_classic() +
  ylab("log2FC") +
  xlab("hours post degron induction") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48))) +
  scale_color_manual(values = colsi,"%GC in CDS")

a<-a + facet_wrap(.~assay,nrow = 1)
a<-a + geom_hline(yintercept = 0,size=.5,lty=2)
pdf(file = "Figures/Fig3E.pdf",width = 7,height = 3.5)
a
dev.off()


df<-df[df$assay%in%c("synthesis","degradation"),]
colsi<-alpha(c("blue","dark red","dark gray"),c(.8,.8,.8))
names(colsi)<-c("TE_down","TE_up","unchanging")

a<-ggplot(df,aes(x=time,y=fcs,color=tx_type_RiboRNA,group=tx_type_RiboRNA))  + stat_summary(geom="line") +stat_summary()

a<-a + theme_classic() +
  ylab("log2FC") +
  xlab("hours post degron induction") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"Tx_class") + scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48)))
a<-a + facet_wrap(.~assay,nrow = 1)
a<-a + geom_hline(yintercept = c(-.5,0,.5),size=.5,lty=2)

pdf(file = "Figures/Fig2A.pdf",width = 8,height = 8)
a + coord_cartesian(ylim = c(-0.5,.5)) + facet_wrap(.~assay,ncol = 1)
dev.off()


dfabc<-df
dfabc<-dfabc[dfabc$assay!="degradation",]
a<-ggplot(dfabc,aes(x=time,y=fcs,color=gcqnt,group=gcqnt))  + stat_summary(geom="line") +stat_summary()

colsi<-alpha(colorRampPalette(colors = c("forestgreen","grey","purple"))(5),c(.8,.8,.8,.8,.8))


colsi<-alpha(c("blue","dark red","dark gray"),c(.8,.8,.8))
names(colsi)<-c("TE_down","TE_up","unchanging")

a<-ggplot(df,aes(x=time,y=fcs,color=gcqnt,group=gcqnt))  + stat_summary(geom="line") +stat_summary()
colsi<-alpha(colorRampPalette(colors = c("forestgreen","grey","purple"))(5),c(.8,.8,.8,.8,.8))

a<-a + theme_classic() +
  ylab("log2FC") +
  xlab("hours post degron induction") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"%GC in CDS")+
  scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48)))
a<-a + facet_wrap(.~assay,nrow = 1)
a<-a + geom_hline(yintercept = 0,size=.5,lty=2)

pdf(file = "Figures/Supplementary_Figure_5.pdf",width = 12,height = 3.5)
a + coord_cartesian(ylim = c(-0.5,.5))
dev.off()


dfaspl<-split(dfa,dfa$experiment)
dfalink<-dfaspl$`04h`
dfalink$posx<-dfalink$RNA_log2FC/2
dfalink$posy<-dfalink$Ribo_log2FC/2
dfalink$dirx<-dfalink$RNA_log2FC/2
dfalink$diry<-dfalink$Ribo_log2FC/2
dfalink$tocol<-dfalink$GCpct_cds*100
#dfalink<-dfalink[dfalink$gene_name%in%c("RHBDD2","ODC1"),]

dfalink<-dfalink[,c("gene_id","posx","posy","dirx","diry","tocol","experiment")]

for(tm in 2:length(dfaspl)){
  dff<-dfaspl[[tm]]
  dff$posx<-dfaspl[[tm-1]]$RNA_log2FC
  dff$posy<-dfaspl[[tm-1]]$Ribo_log2FC
  dff$tocol<-dfaspl[[tm]]$GCpct_cds*100
  
  dff$dirx<-dff$RNA_log2FC-dff$posx
  dff$diry<-dff$Ribo_log2FC-dff$posy
  #dff<-dff[dff$gene_name%in%c("RHBDD2","ODC1"),]
  dff<-dff[,c("gene_id","posx","posy","dirx","diry","tocol","experiment")]
  dfalink<-rbind(dfalink,dff)
}
gnsok<-dfa[which(dfa$RiboRNA_padj<.01 | dfa$RNA_padj<.01 | dfa$Ribo_padj<.01),"gene_id"]
dfalink<-dfalink[dfalink$gene_id%in%gnsok,]
resol=70
axm<-max(abs(c(dfa$RNA_log2FC,dfa$Ribo_log2FC)),na.rm = T)
bns<-seq(-axm,axm,length.out = resol+1)
grida<-cbind(bns[-length(bns)],bns[-1])
gridp<-paste(grida[,1],grida[,2],sep = ";")
gridp<-expand.grid(gridp,gridp)
gridp<-paste(gridp[,1],gridp[,2],sep = ";")
grida<-t(sapply(strsplit(gridp,";"),c))
mode(grida)<-"numeric"

#resiall$tocol<-scale(resiall$tocol)

resiall<-c()
for(i in unique(dfalink$experiment)){
    resi_x<-matrix(0,nrow = resol,ncol = resol)
    resi_y<-resi_x
    resi_n<-resi_x
    resi_tocol<-resi_x
    
    dff<-dfalink[dfalink$experiment==i,]
    #gridx<-apply(dff,1,function(x){which(grida[,1]<=x["posx"] & grida[,2]>=x["posx"] & grida[,3]<=x["posy"] &grida[,4]>=x["posy"])})
    
    posx<-as.numeric(cut(dff$posx,breaks = bns,include.lowest = T,right = T))
    posy<-as.numeric(cut(dff$posy,breaks = bns,include.lowest = T,right = T))
    idx<-paste(posx,posy,sep = ";")
    posssx<-NumericList(split(dff$dirx,idx))
    coor<-t(sapply(strsplit(names(posssx),";"),c))
    mode(coor)<-"integer"
    posssy<-mean(NumericList(split(dff$diry,idx)),na.rm=T)
    tocol<-mean(NumericList(split(dff$tocol,idx)),na.rm=T)
    
    resi_n[coor]<-elementNROWS(posssx)
    resi_x[coor]<-unname(mean(posssx))
    resi_y[coor]<-unname(posssy)
    resi_tocol[coor]<-unname(tocol)
    ressi<-cbind.data.frame(c(resi_n),c(resi_x),c(resi_y),c(resi_tocol))
    colnames(ressi)<-c("obs","dirx","diry","tocol")
    ressi$experiment<-i
    
    resiall<-rbind(resiall,ressi)
}

resiall$posx<-rep(grida[,1]+diff(bns)[1]/2,length(unique(dfalink$experiment)))
resiall$posy<-rep(grida[,3]+diff(bns)[1]/2,length(unique(dfalink$experiment)))
resiall_cp<-resiall

#resiall<-resiall[resiall$obs>0,]
#resiall$obs[resiall$obs>50]<-50

#resiall<-resiall[resiall$experiment%in%c("04h_vs_DMSO","08h_vs_DMSO","16h_vs_DMSO"),]
#resiall$tocol<-scale(resiall$tocol)
#resiall$tocol[resiall$tocol<(-.5)]<-(-.5)
#resiall$tocol[resiall$tocol>(.5)]<-.5
invix<-function(x){1-1/x}
xivni<-function(x){-(1/(x-1))}

resiall<-resiall_cp
resiall$obs2<-resiall$obs
#resiall$obs2[resiall$obs>20]<-20
resiall<-resiall[resiall$obs>2,]
a<-ggplot(resiall,aes(x=posx,y=posy))
a<-a + theme_bw() +
    geom_segment(aes(x=posx,xend=posx+dirx,y=posy,yend=posy+diry,alpha=obs2,color=tocol),size=1,arrow = arrow(angle = 20, length = unit(.2, "cm"), type = "open")) +
    #geom_streamline(aes(x=posx,dx=dirx,y=posy,dy=diry),res=10) +
    scale_color_gradient2(low = "forestgreen",mid = "grey",high = "purple",midpoint = mean(resiall$tocol),"%GC in CDS")+
    #scale_color_gradient2(low = "dark blue",mid = "light grey",high = "dark orange",midpoint = 0)+
    ylab(paste("Ribo-seq log2FC")) +
    xlab(paste("RNA-seq log2FC")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
    #scale_color_manual(values = colsi,"Tx_class") +
    #scale_alpha_continuous("# mRNAs",trans = scales::trans_new(name = "test",transform =invix,inverse = xivni))
    scale_alpha_continuous("# mRNAs",trans = scales::trans_new(name = "test",transform =invix,inverse = xivni),range = c(0.2,1),breaks=c(5,10,50))
a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
#ariborna<-a + ylim(-axm,axm) + xlim(-axm,axm)
#pdf(file = "Figures//bfd/lcalviel/data/riboseq/DDX3_degrtime/DEDEX/vector_pl_ht.pdf",width = 20,height = 5)
axm=2
a<-a + coord_cartesian(xlim = c(-axm,axm),ylim=c(-axm,axm)) + facet_wrap(experiment~.,scales = "fixed",nrow =  1)
#dev.off()
pdf(file = "Figures/Fig3D.pdf",width = 22,height = 5)
a
dev.off()

# make a vector plot for Fig 1

dfaspl<-split(dfa,dfa$experiment)
dfalink<-dfaspl[["04h"]]
dfalink$posx<-dfalink$RNA_log2FC/2
dfalink$posy<-dfalink$Ribo_log2FC/2
dfalink$dirx<-dfalink$RNA_log2FC/2
dfalink$diry<-dfalink$Ribo_log2FC/2
dfalink$tocol<-dfaspl[["04h"]]$delta_TE
#dfalink<-dfalink[dfalink$gene_name%in%c("RHBDD2","ODC1"),]

dfalink<-dfalink[,c("gene_id","posx","posy","dirx","diry","tocol","experiment")]

for(tm in 2:length(dfaspl)){
  dff<-dfaspl[[tm]]
  dff$posx<-dfaspl[[tm-1]]$RNA_log2FC
  dff$posy<-dfaspl[[tm-1]]$Ribo_log2F
  dff$tocol<-dfaspl[[tm-1]]$delta_TE
  #dff$tocol<-dfaspl$time48h$delta_TE
  
  dff$dirx<-dff$RNA_log2FC-dff$posx
  dff$diry<-dff$Ribo_log2FC-dff$posy
  #dff<-dff[dff$gene_name%in%c("RHBDD2","ODC1"),]
  dff<-dff[,c("gene_id","posx","posy","dirx","diry","tocol","experiment")]
  dfalink<-rbind(dfalink,dff)
}
gnsok<-dfa[which(dfa$RiboRNA_padj<.01 | dfa$RNA_padj<.01 | dfa$Ribo_padj<.01),"gene_id"]
dfalink<-dfalink[dfalink$gene_id%in%gnsok,]
resol=70
axm<-max(abs(c(dfa$RNA_log2FC,dfa$Ribo_log2FC)),na.rm = T)
bns<-seq(-axm,axm,length.out = resol+1)
grida<-cbind(bns[-length(bns)],bns[-1])
gridp<-paste(grida[,1],grida[,2],sep = ";")
gridp<-expand.grid(gridp,gridp)
gridp<-paste(gridp[,1],gridp[,2],sep = ";")
grida<-t(sapply(strsplit(gridp,";"),c))
mode(grida)<-"numeric"


resiall<-c()
for(i in unique(dfalink$experiment)){
  resi_x<-matrix(0,nrow = resol,ncol = resol)
  resi_y<-resi_x
  resi_n<-resi_x
  resi_tocol<-resi_x
  
  dff<-dfalink[dfalink$experiment==i,]
  #gridx<-apply(dff,1,function(x){which(grida[,1]<=x["posx"] & grida[,2]>=x["posx"] & grida[,3]<=x["posy"] &grida[,4]>=x["posy"])})
  
  posx<-as.numeric(cut(dff$posx,breaks = bns,include.lowest = T,right = T))
  posy<-as.numeric(cut(dff$posy,breaks = bns,include.lowest = T,right = T))
  idx<-paste(posx,posy,sep = ";")
  posssx<-NumericList(split(dff$dirx,idx))
  coor<-t(sapply(strsplit(names(posssx),";"),c))
  mode(coor)<-"integer"
  posssy<-mean(NumericList(split(dff$diry,idx)),na.rm=T)
  tocol<-mean(NumericList(split(dff$tocol,idx)),na.rm=T)
  
  resi_n[coor]<-elementNROWS(posssx)
  resi_x[coor]<-unname(mean(posssx))
  resi_y[coor]<-unname(posssy)
  resi_tocol[coor]<-unname(tocol)
  ressi<-cbind.data.frame(c(resi_n),c(resi_x),c(resi_y),c(resi_tocol))
  colnames(ressi)<-c("obs","dirx","diry","tocol")
  ressi$experiment<-i
  
  resiall<-rbind(resiall,ressi)
}

resiall$posx<-rep(grida[,1]+diff(bns)[1]/2,length(unique(dfalink$experiment)))
resiall$posy<-rep(grida[,3]+diff(bns)[1]/2,length(unique(dfalink$experiment)))
resiall_cp<-resiall

#resiall<-resiall[resiall$obs>0,]
#resiall$obs[resiall$obs>50]<-50

#resiall<-resiall[resiall$experiment%in%c("04h_vs_DMSO","08h_vs_DMSO","16h_vs_DMSO"),]
#resiall$tocol<-scale(resiall$tocol)
#resiall$tocol[resiall$tocol<(-.5)]<-(-.5)
#resiall$tocol[resiall$tocol>(.5)]<-.5
invix<-function(x){1-1/x}
xivni<-function(x){-(1/(x-1))}

resiall<-resiall_cp
resiall$obs2<-resiall$obs
#resiall$obs2[resiall$obs>20]<-20
resiall<-resiall[resiall$obs>2,]



a<-ggplot(resiall,aes(x=posx,y=posy))
a<-a + theme_bw() +
  geom_segment(aes(x=posx,xend=posx+dirx,y=posy,yend=posy+diry,alpha=obs2,color=tocol),size=1,arrow = arrow(angle = 20, length = unit(.2, "cm"), type = "open")) +
  #geom_streamline(aes(x=posx,dx=dirx,y=posy,dy=diry),res=10) +
  scale_color_gradient2(low = "dark blue",mid = "grey",high = "dark red",midpoint = 0,"delta_TE")+
  #scale_color_gradient2(low = "dark blue",mid = "light grey",high = "dark orange",midpoint = 0)+
  ylab(paste("Ribo-seq log2FC")) +
  xlab(paste("RNA-seq log2FC")) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
  #scale_color_manual(values = colsi,"Tx_class") +
  #scale_alpha_continuous("# mRNAs",trans = scales::trans_new(name = "test",transform =invix,inverse = xivni))
  scale_alpha_continuous("# mRNAs",trans = scales::trans_new(name = "test",transform =invix,inverse = xivni),range = c(0.2,1),breaks=c(5,10,50))
a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
#ariborna<-a + ylim(-axm,axm) + xlim(-axm,axm)
#pdf(file = "Figures//bfd/lcalviel/data/riboseq/DDX3_degrtime/DEDEX/vector_pl_ht.pdf",width = 20,height = 5)
axm=2
a<-a + coord_cartesian(xlim = c(-axm,axm),ylim=c(-axm,axm)) + facet_wrap(experiment~.,scales = "fixed",nrow =  1)
#dev.off()
pdf(file = "Figures/Fig1D.pdf",width = 22,height = 5)
a
dev.off()



dfa<-scatterplot_RiboRNAIntronExon(summary_multiDE)[[2]]$data
dfa<-dfa[!is.na(dfa$xy_RiboRNA),]
dfa<-dfa[dfa$experiment=="48h_vs_DMSO",]
dfa$xy_RiboRNA<-as.character(dfa$xy_RiboRNA)
dfa$xy_RiboRNA[dfa$xy_RiboRNA%in%c("xy","x")]<-"other"
#dfa$xy_RiboRNA[which(dfa$xy_RiboRNA%in%c("x","xy"))]<-"other"
dfa$xy_RiboRNA<-factor(dfa$xy_RiboRNA,levels=c("-xy","y","other"))
#only some xy, check tomorrow
collxy<-alpha(c("purple","blue","dark grey"),alpha = .6)
names(collxy)<-levels(dfa$xy_RiboRNA)

a<-ggplot(dfa,aes(x=RNA_log2FC,y=Ribo_log2FC,color=xy_RiboRNA,size=RiboRNA_padj)) + geom_point()

a<-a + theme_bw() +
  ylab(paste("Ribo-seq "," log2FC",sep = "")) +
  xlab(paste("RNA-seq log2FC",sep = "")) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values=collxy,"") +
  scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv,
                                                                               domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))

axm=3.5
a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + coord_cartesian(xlim = c(-axm,axm),ylim=c(-axm,axm))
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1)
pdf(file = "Figures/Fig3A.pdf",width = 10,height = 6)
ariborna
dev.off()


#load("data/degrtime_multiDE_results_RF.RData")
load("data/degron_multiDE_results_RFnew.RData")

#rfcor<-get_rf_corrs(res_dfs)
#rfimpor<-get_rf_importance(res_dfs)
#rflscomp<-compare_lasso_rf_importance(res_dfs)


list_corplts<-list()
rfone<-res_dfs
nims<-names(rfone)
nims<-nims[nims!="df_list"]
for(j in nims){
  rfone2<-rfone[[j]]
  nims2<-names(rfone2)
  ress_rf2<-list()
  for(i in nims2){
    res_corrs<-lapply(rfone2[[i]],function(x){
      xx<-t(sapply(x,function(y){as.numeric(y$res_tests[1:4])}))
      colnames(xx)<-names(x[[1]]$res_tests[1:4])
      xx
    })
    
    res_corrs_lass<-lapply(rfone2[[i]],function(x){
      xx<-t(sapply(x,function(y){as.numeric(y$res_tests[6:9])}))
      colnames(xx)<-names(x[[1]]$res_tests[1:4])
      xx
    })
    
    ress<-data.frame(do.call(res_corrs,what = rbind))
    ress$feat<-rep(names(res_corrs),each=dim(ress)[1]/length(res_corrs))
    ress<-suppressMessages(melt(ress))
    ress<-ress[grep(ress$feat,pattern = "posdens",invert = T),]
    ress$method="Pearson"
    ress$method[grepl(ress$variable,pattern = "corrs")]="Spearman"
    ress$dataset="train"
    ress$dataset[grepl(ress$variable,pattern = "tst")]="test"
    ress$xy<-i
    
    ress$algo="Random_forest"
    
    
    ress_rf<-ress
    ress<-data.frame(do.call(res_corrs_lass,what = rbind))
    ress$feat<-rep(names(res_corrs),each=dim(ress)[1]/length(res_corrs))
    ress<-suppressMessages(melt(ress))
    ress<-ress[grep(ress$feat,pattern = "posdens",invert = T),]
    ress$method="Pearson"
    ress$method[grepl(ress$variable,pattern = "corrs")]="Spearman"
    ress$dataset="train"
    ress$dataset[grepl(ress$variable,pattern = "tst")]="test"
    ress$xy<-i
    ress$algo="lasso"
    ress<-rbind(ress_rf,ress)
    
    ress_rf2[[i]]<-ress
  }
  ress<-do.call(ress_rf2,what = rbind)
  ress$algo<-factor(ress$algo,levels=c("Random_forest","lasso"))
  #ress$xy<-factor(ress$xy,levels=c("all","RiboRNA_y","RiboRNA_xy","RiboRNA_-xy","RiboRNA_x","IntronExon_y","IntronExon_xy","IntronExon_-xy","IntronExon_x"))
  
  #ress_rfimp<-ress_rfimp[grep(ress_rfimp$xy,pattern = "xy",invert = T),]
  
  list_corplts[[j]]<-ress
}

ress<-list_corplts$`48h`
ress<-ress[ress$xy%in%c("all","RiboRNA_y","RiboRNA_-xy","RiboRNA_x","RiboRNA_xy"),]
ress<-ress[ress$feat=="delta_TE",]
ress$xy<-factor(ress$xy,levels=c("all","RiboRNA_y","RiboRNA_-xy","RiboRNA_x","RiboRNA_xy"))


ress$xy<-gsub(ress$xy,pattern = "RiboRNA_",replacement = "")
ress$xy<-factor(ress$xy,levels=c("all","y","-xy","x","xy"))

collxy<-alpha(c("black","blue","purple","mediumorchid4","dark grey"),alpha = .6)
names(collxy)<-levels(ress$xy)

ress<-ress[ress$xy%in%c("all","y","-xy"),]
rf_corplot<-ggplot(ress,aes(y=value,x=feat,shape=dataset,color=xy)) +
  #geom_jitter() +
  stat_summary(position=position_dodge(.9),size=.7) +
  ylab("Correlation\nobserved-predicted") +
  xlab("") +
  scale_color_manual(values = collxy,"") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))+
  facet_wrap(method~algo) + coord_cartesian(ylim = c(0,1))+
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.1, color="black" )
  )

#pdf(file = "Figures/Fig2_supp_corr_rfls.pdf",width = 11,height = 6)
#rf_corplot
#dev.off()


collxy<-alpha(c("black","blue","purple","mediumorchid4","dark grey"),alpha = .6)
names(collxy)<-levels(ress$xy)

ress$xy<-factor(as.character(ress$xy),levels=c("all","y","-xy"))
collxy<-collxy[levels(ress$xy)]
ress<-ress[ress$variable=="corrp_tst" & ress$algo=="Random_forest",]
rf_corplot<-ggplot(ress,aes(y=value,x=xy,color=xy)) +
  #geom_jitter() +
  stat_summary(position=position_dodge(.9),size=.7) +geom_jitter(position=position_dodge(.9),size=.4)+
  ylab("Correlation with test data\npredicted vs. observed") +
  xlab("") +
  scale_color_manual(values = collxy,"") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))+
  coord_cartesian(ylim = c(0,1))+
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.1, color="black" )
  )


pdf(file = "Figures/Fig3B.pdf",width = 6,height = 5)
rf_corplot
dev.off()

nims<-names(rfone)
nims<-nims[nims!="df_list"]
j="48h"
rfone2<-rfone[[j]]
nims2<-names(rfone2)
ress_rfimp<-list()
for(i in nims2){
  
  ress<-lapply(rfone2[[i]],function(x){
    xx<-melt(lapply(x,function(y){y$importance}))
    xx<-xx[xx$Var2=="%IncMSE",]
    xx$value<-xx$value/sum(xx$value)
    xx$xy<-i
    xx
  })
  
  for(aa in names(ress)){
    ress[[aa]]$feature<-aa
  }
  ress_imp<-do.call(ress,what = rbind)
  ress_rfimp[[i]]<-ress_imp
}
ress_rfimp<-do.call(ress_rfimp,what = rbind)
ress_rfimp<-ress_rfimp[ress_rfimp$xy%in%c("all","RiboRNA_y","RiboRNA_-xy"),]

ress_rfimp<-ress_rfimp[ress_rfimp$feature=="delta_TE",]

topsi<-rev(as.character(unique(ress_rfimp$Var1[order(ress_rfimp$value,decreasing = T)])))

#ress_rfimp$Var1<-factor(ress_rfimp$Var1,levels=unique(ress_rfimp$Var1))
ress_rfimp$Var1<-factor(ress_rfimp$Var1,levels=topsi)

vars<-unique(ress_rfimp$Var1)
ress_rfimp$xy<-factor(ress_rfimp$xy,levels=c("all","RiboRNA_y","RiboRNA_xy","RiboRNA_-xy","RiboRNA_x","IntronExon_y","IntronExon_xy","IntronExon_-xy","IntronExon_x"))

collxy<-alpha(c("black","blue","steelblue","purple","dark blue","red","darkolivegreen2","darkolivegreen4","dark red"),alpha = .6)
names(collxy)<-levels(ress_rfimp$xy)
#ress_rfimp<-ress_rfimp[grep(ress_rfimp$xy,pattern = "xy",invert = T),]
set.seed(667)
clusti<-kmeans(ress_imp$value,centers = 3)$cluster


cl_high<-unname(clusti[which.max(ress_imp$value)])
cl_low<-unname(clusti[which.min(ress_imp$value)])
cl_mid<-unique(clusti)[!unique(clusti)%in%c(cl_high,cl_low)]
okvarse<-unique(c(as.character(ress_imp$Var1[clusti==cl_high]),as.character(ress_imp$Var1[clusti==cl_mid])))
aaa<-ifelse(!unique(ress_imp$Var1)%in%okvarse, "grey", "black")

#levels(ress_rfimp$Var1)<-as.character(unique(ress_rfimp$Var1[order(ress_rfimp$value,decreasing = T)]))

rf_impor<-ggplot(ress_rfimp,aes(y=Var1,x=value,color=xy)) + geom_jitter(size=.2)+ stat_summary() +theme_bw()+
  ylab("") +
  scale_color_manual(values = collxy,"") +
  xlab("scaled feature Importance") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  #theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=7,colour=aaa))  +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=7))  +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  theme(strip.text.x = element_text(size = 15,face = "bold"))
rf_impor<-rf_impor+facet_wrap(~feature,nrow=1,scale="free_x")


new_ress_rfimp<-ress_rfimp
rf_others<-levels(new_ress_rfimp$Var1)[1:85]
new_ress_rfimp$Var1<-as.character(new_ress_rfimp$Var1)
new_ress_rfimp$Var1[new_ress_rfimp$Var1%in%rf_others]<-"other"
new_ress_rfimp<-aggregate(new_ress_rfimp$value,list(new_ress_rfimp$Var1,new_ress_rfimp$L1,new_ress_rfimp$xy,new_ress_rfimp$feature),mean)
colnames(new_ress_rfimp)<-c("Var1","L1","xy","feature","value")
topsi<-rev(as.character(unique(new_ress_rfimp$Var1[order(new_ress_rfimp$value,decreasing = T)])))
#ress_rfimp$Var1<-factor(ress_rfimp$Var1,levels=unique(ress_rfimp$Var1))
new_ress_rfimp$Var1<-factor(new_ress_rfimp$Var1,levels=topsi)

new_ress_rfimp$xy<-as.character(new_ress_rfimp$xy)
new_ress_rfimp$xy<-gsub(new_ress_rfimp$xy,pattern = "RiboRNA_",replacement = "")
new_ress_rfimp$xy<-factor(new_ress_rfimp$xy,levels=c("all","y","-xy"))


collxy<-alpha(c("black","blue","purple","mediumorchid4","dark grey"),alpha = .6)
names(collxy)<-levels(new_ress_rfimp$xy)
collxy<-collxy[levels(new_ress_rfimp$xy)]

rf_impor<-ggplot(new_ress_rfimp,aes(y=Var1,x=value,color=xy)) + geom_jitter(size=.2)+ stat_summary() +theme_bw()+
  ylab("") +
  scale_color_manual(values = collxy,"") +
  xlab("scaled feature Importance") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  #theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=7,colour=aaa))  +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=7))  +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  theme(strip.text.x = element_text(size = 15,face = "bold"))
rf_impor<-rf_impor+facet_wrap(~feature,nrow=1,scale="free_x")

pdf(file = "Figures/Fig3C.pdf",width = 6,height = 5)
rf_impor
dev.off()



nims<-names(rfone)
nims<-nims[nims!="df_list"]

rfone2<-rfone[[j]]
nims2<-names(rfone2)
ress_rfimp<-list()
for(i in nims2){
  
  ress<-lapply(rfone2[[i]],function(x){
    xx<-melt(lapply(x,function(y){y$importance}))
    xx<-xx[xx$Var2=="%IncMSE",]
    xx$value<-xx$value/sum(xx$value)
    xx$xy<-i
    xx
  })
  
  for(aa in names(ress)){
    ress[[aa]]$feature<-aa
  }
  ress_imp<-do.call(ress,what = rbind)
  #ress_imp$Var1<-NULL
  ress_imp$Var2<-NULL
  
  ress_imp_rf<-ress_imp
  ress<-lapply(rfone2[[i]],function(x){
    xx<-melt(lapply(x,function(y){y$lasso_coefs}))
    xx$Var1<-rep(names(x[[1]]$lasso_coefs),dim(xx)[1]/length(names(x[[1]]$lasso_coefs)))
    xx<-xx[,c(3,1,2)]
    #xx$value<-xx$value/sum(xx$value)
    xx<-xx[xx$Var1!="(Intercept)",]
    xx$xy<-i
    xx
  })
  
  for(aa in names(ress)){
    ress[[aa]]$feature<-aa
  }
  ress_imp<-do.call(ress,what = rbind)
  
  ress_imp$value_rf<-ress_imp_rf$value
  ress_rfimp[[i]]<-ress_imp
}
ress_rfimp<-do.call(ress_rfimp,what = rbind)

ress_rfimp$Var1<-factor(ress_rfimp$Var1,levels=unique(ress_rfimp$Var1))

vars<-unique(ress_rfimp$Var1)

ress_rfimp$xy<-factor(ress_rfimp$xy,levels=c("all","RiboRNA_y","RiboRNA_-xy"))

collxy<-alpha(c("black","blue","purple","steelblue","dark blue","red","darkolivegreen2","darkolivegreen4","dark red"),alpha = .6)
names(collxy)<-levels(ress_rfimp$xy)
ress_rfimp<-ress_rfimp[ress_rfimp$xy%in%c("all","RiboRNA_y","RiboRNA_-xy"),]
ress_rfimp<-ress_rfimp[ress_rfimp$feature=="delta_TE",]

#minvalo<-apply(cbind(abs(ress_rfimp$value),ress_rfimp$value_rf),1,mean)
minvalo<-ress_rfimp$value_rf
set.seed(666)
clusti<-kmeans(minvalo,centers = 3)$cluster


cl_high<-unname(clusti[which.max(minvalo)])
cl_low<-unname(clusti[which.min(minvalo)])
cl_mid<-unique(clusti)[!unique(clusti)%in%c(cl_high,cl_low)]
okvarse<-unique(c(as.character(ress_rfimp$Var1[clusti==cl_high]),as.character(ress_rfimp$Var1[clusti==cl_mid])))
if(length(okvarse)>10){okvarse=unique(ress_rfimp$Var1[order(ress_rfimp$value_rf,decreasing = T)])[1:10]}


mns<-aggregate(ress_rfimp$value,by=list(ress_rfimp$Var1,ress_rfimp$feature,ress_rfimp$xy),function(x){c(mean(x),sd(x))})
mns$mean_L<-mns$x[,1]
mns$sd_L<-mns$x[,2]
mns$x<-NULL
mns<-mns[order(paste(mns$Group.1,mns$Group.2,mns$Group.3)),]

aaa2<-mns$mean_L>0

mns<-aggregate(abs(ress_rfimp$value),by=list(ress_rfimp$Var1,ress_rfimp$feature,ress_rfimp$xy),function(x){c(mean(x),sd(x))})
mns$mean_L<-mns$x[,1]
mns$sd_L<-mns$x[,2]
mns$x<-NULL
mns<-mns[order(paste(mns$Group.1,mns$Group.2,mns$Group.3)),]

all_sdm<-mns
mns<-aggregate(ress_rfimp$value_rf,by=list(ress_rfimp$Var1,ress_rfimp$feature,ress_rfimp$xy),function(x){c(mean(x),sd(x))})
mns$mean_RF<-mns$x[,1]
mns$sd_RF<-mns$x[,2]
mns$x<-NULL
mns<-mns[order(paste(mns$Group.1,mns$Group.2,mns$Group.3)),]

all_sdm$mean_rf<-mns$mean_RF
all_sdm$sd_rf<-mns$sd_RF

vars<-unique(all_sdm$Group.1)
all_sdm$labbe<-NA
for(ixo in unique(all_sdm$Group.2)){
  
  for(ixo2 in all_sdm$Group.3){
    uici<-all_sdm$Group.2==ixo & all_sdm$Group.3==ixo2
    vare<-unique(c(as.character(all_sdm$Group.1[uici][order(all_sdm$mean_rf[uici],decreasing = T)][1:5]),as.character(all_sdm$Group.1[uici][order(all_sdm$mean_L[uici],decreasing = T)][1:5])))
    all_sdm$labbe[all_sdm$Group.1%in%vare & uici]<-as.character(all_sdm$Group.1[all_sdm$Group.1%in%vare & uici])
  }
}
#all_sdm$labbe[!all_sdm$labbe%in%okvarse]<-NA

all_impor<-ggplot(all_sdm,aes(y=mean_L,x=mean_rf,color=aaa2,label=labbe)) +  geom_point() +theme_bw()+
  ylab("Feature importance (Lasso)\nabsolute values") +
  xlab("Feature Importance (Random Forest)") +
  scale_color_manual(values = c("red","blue"),"Positive\npredictor") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  #scale_color_manual(values = aaa2)+
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  theme(strip.text.x = element_text(size = 15,face = "bold")) + theme(legend.position="none")+
  geom_errorbar(aes(ymin = mean_L-sd_L,ymax = mean_L+sd_L)) +
  geom_errorbarh(aes(xmin = mean_rf-sd_rf,xmax = mean_rf+sd_rf)) +
  geom_text_repel(size=5,force=32,show.legend = F)
all_impor<-all_impor+facet_wrap(.~Group.3,nrow = 1,scales = "free_x")


pdf(file = "Figures/Supplementary_Figure_3.pdf",width = 20,height = 7)
all_impor
dev.off()


aa<-rfone$df_list$`48h`
aa<-aa[which(aa$TPM_RNA>3 & !is.na(aa$delta_TE)),]
aa$GCpct_cds<-aa$GCpct_cds*100
aa<-aa[!is.infinite(aa$delta_TE),]
a<-ggplot(aa,aes(x=delta_TE,y=GCpct_cds)) + geom_point()
a<-a + theme_bw() +
  ylab("%GC in CDS") +
  xlab("delta TE") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))

#pdf(file = "Figures/Fig2_time48_supp_gc.pdf",width = 10,height = 6)
#a + geom_smooth(method="lm")
#dev.off()    

aa<-rfone$df_list$`48h`
aa<-aa[which(aa$TPM_RNA>3 & !is.na(aa$delta_TE)),]
aa$GCpct_cds<-aa$GCpct_cds*100
aa<-aa[!is.infinite(aa$delta_TE),]
a<-ggplot(aa,aes(x=RNA_log2FC,y=GCpct_cds)) + geom_point()
a<-a + theme_bw() +
  ylab("%GC in CDS") +
  xlab("RNA log2FC") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))

#pdf(file = "Figures/Fig2_time48_supp_gc2.pdf",width = 10,height = 6)
#a + geom_smooth(method="lm")
#dev.off()                                                                    


load("data/new_rf_resoverview")
load("data/new_multiDE_resoverview")
#rfres<-do.call(list_rfres,what=rbind)
ginon<-GTF_annotation$trann$gene_name[match(rownames(enc_dfres[[1]]),GTF_annotation$trann$gene_id)]

#arbipi<-sapply(strsplit(rfres$condition,"_"),"[[",1)
#cell<-sapply(strsplit(rfres$condition,"_"),"[[",2)
#meth<-sapply(strsplit(rfres$condition,"_"),"[[",3)

ginon<-summary_multiDE[[1]][[1]]$gene_name

arbipi<-sapply(strsplit(names(enc_dfres),"_"),"[[",1)
cell<-sapply(strsplit(names(enc_dfres),"_"),"[[",2)
meth<-sapply(strsplit(names(enc_dfres),"_"),"[[",3)

kds<-c()
pvs<-c()
n_sigi<-c()
for(i in 1:length(enc_dfres)){
  aa<-enc_dfres[[i]]
  kds<-c(kds,aa$RNA_log2FC[ginon==arbipi[i]])
  pvs<-c(pvs,aa$RNA_padj[ginon==arbipi[i]])
  n_sigi<-c(n_sigi,length(which(aa$RNA_padj<.01)))
}

kdos<-kds

#G3BP2 very different bwt cell lines
names(kdos)<-names(enc_dfres)
okko<-n_sigi>100 & kds<(-1)
kdos<-kdos[okko]

#kdos<-kdos[grep(names(kdos),pattern = "K562")]

best_kds<-unname(sapply(split(kdos,paste(arbipi,cell,meth,sep = "_")[okko]),function(x){names(x)[which.min(x)]}))
kdos<-kdos[best_kds]

enc_rfres$RBP<-sapply(strsplit(enc_rfres$condition,"_"),"[[",1)
enc_rfres$cell<-sapply(strsplit(enc_rfres$condition,"_"),"[[",2)
enc_rfres$method<-sapply(strsplit(enc_rfres$condition,"_"),"[[",3)

dfcorrs<-enc_rfres[enc_rfres$condition%in%best_kds,]

#dfcorrs<-rfres
dfcorrs<-unique(dfcorrs[,c("RBP","cell","method","xy","feature","corrp_tst")])
dfcorrs<-dfcorrs[dfcorrs$feature=="RNA_log2FC",]
dfcorrs_cp<-dfcorrs

dfcorrs<-dfcorrs[dfcorrs$xy=="all",]

dfcorrs<-dfcorrs[dfcorrs$method=="shRNA",]
dfcorrs<-dfcorrs[dfcorrs$cell=="K562",]
dfcorrs<-dfcorrs[order(dfcorrs$corrp_tst,decreasing = T),]
dfcorrs$type_exp<-interaction(dfcorrs$cell,dfcorrs$method)
ago<-aggregate(dfcorrs$corrp_tst,by=list(dfcorrs$RBP),mean)
dfcorrs_cp<-dfcorrs
aribipisi<-ago[,1][order(ago[,2],decreasing = T)]
dfcorrs$RBP<-factor(dfcorrs$RBP,levels = aribipisi)
colsi<-setNames(object = alpha(c("black","blue","purple","red","orange"),.8),nm = c("all","x","-xy","y","xy"))
a<-ggplot(dfcorrs,aes(x=RBP,y=corrp_tst)) + stat_summary()
a<-a + theme_classic() +
  xlab("") +
  ylab("correlation with test data") +
  theme(axis.title.x = element_text(size=12),axis.text.x  = element_text(angle=45, vjust=0.5, size=10)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  
#scale_color_manual(values = colsi,"IntronExon_coords")
#pdf(file = "Figures/Fig3_supp_ENC_rf_overview.pdf",width = 20,height = 5)
#a+facet_wrap(feature~.,nrow = 2)
#a
#dev.off()


df2<-dfcorrs_cp
#df2<-df2[df2$method=="shRNA",]
#df2<-df2[df2$cell=="K562",]
df2<-df2[order(df2$corrp_tst,decreasing = T),]
df2$type_exp<-interaction(df2$cell,df2$method)
df2<-df2[df2$feature=="RNA_log2FC",]
df2$condition<-sapply(strsplit(rownames(df2),"[.]"),"[[",1)
df2<-df2[which(df2$xy%in%c("all","x","y")),]
#df2<-df2[df2$RBP%in%c("DDX3X","G3BP1","G3BP2","CNOT10"),]
df2$RBP<-factor(df2$RBP,levels = unique(df2$RBP[df2$xy=="all"]))
df2$xy[df2$xy=="x"]<-"mRNA"
df2$xy[df2$xy=="y"]<-"premRNA"
colsi<-setNames(object = alpha(c("black","blue","purple","red","orange"),.8),nm = unique(df2$xy))
a<-ggplot(df2,aes(x=RBP,y=corrp_tst,color=xy)) + stat_summary()
a<-a + theme_classic() +
  xlab("") +
  ylab("Correlation with test data") +
  theme(axis.title.x = element_text(size=12),axis.text.x  = element_text(angle=45, vjust=0.5, size=10)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"predominant change") + geom_hline(yintercept = c(.25,.5,.75),linetype=2)

#pdf(file = "Figures/Fig3_supp_ENC_rf_overview2.pdf",width = 50,height = 9)
#a+facet_wrap(feature~.,nrow = 2)
#a+facet_wrap(type_exp~.,nrow = length(unique(df2$type_exp)),drop = T)
#dev.off()



dfd<-enc_rfres[enc_rfres$condition%in%best_kds,]

dfd<-dfd[dfd$feature=="RNA_log2FC",]
dfd<-dfd[grep(dfd$Var1,pattern = "ibo|base_TE",invert = T),]
dfd$Var1<-as.character(dfd$Var1)
dfd<-aggregate(list(dfd$value,dfd$corrp_tst),by=list(dfd$Var1,dfd$RBP,dfd$cell,dfd$method,dfd$xy),mean)

colnames(dfd)<-c("Var1","RBP","cell","method","xy","importance","corrp_tst")
dfd<-dfd[dfd$xy=="all",]
corsi<-c()
for(i in unique(dfd$Var1)){
  corsi[i]<-cor.test(dfd$importance[dfd$Var1==i],dfd$corrp_tst[dfd$Var1==i])$estimate
}
df<-dfd[dfd$Var1%in%c("GCpct_cds","GCpct_threeut","codonfr_GCG"),]
df$label<-NA
oka<-df$corrp_tst>.7 & df$importance>.06
df$label[oka]<-df$RBP[oka]
colsi<-setNames(object = alpha(c("black","blue","purple","red","orange"),.8),nm = unique(df2$xy))
df$experiment<-interaction(df$cell,df$method)
a<-ggplot(df,aes(x=importance,y=corrp_tst,color=Var1,label=label,shape=experiment)) + geom_smooth() + geom_point(size=.9)
a<-a + theme_classic() +
  xlab("Feature importance") +
  ylab("Correlation with test data") +
  theme(axis.title.x = element_text(size=24),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = c("dark grey","dark blue","steelblue"),"Feature")

pdf(file = "Figures/Fig4C.pdf",width = 10,height = 6)
a + geom_text_repel(show.legend = F)
dev.off()

aa<-enc_dfres$DDX3X_K562_shRNA_version1
#aa<-data.frame(cbind(enc_dfres$DDX3X_HepG2_shRNA_version1$RNA_log2FC[enc_dfres$DDX3X_HepG2_shRNA_version1$TPM_RNA>3],summary_multiDE$uniq[[1]]$GCpct_cds[enc_dfres$DDX3X_HepG2_shRNA_version1$TPM_RNA>3]))
aa$GCpct_cds<-summary_multiDE$uniq[[1]]$GCpct_cds
aa<-aa[complete.cases(aa),]

#colnames(aa)<-c("RNA_log2FC","GCpct_cds")
aa$GCpct_cds<-aa$GCpct_cds*100
a<-ggplot(aa,aes(x=RNA_log2FC,y=GCpct_cds,size=RNA_padj)) + geom_point(alpha=.5)
a<-a + theme_bw() +
  ylab("%GC in CDS") +
  xlab("RNA log2FC") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))+
  scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv,
                                                                 domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))

pdf(file = "Figures/Fig4D.pdf",width = 10,height = 6)
a + geom_smooth(show.legend = F) + coord_cartesian(xlim = c(-2.25,2.25),ylim=c(25,85))
dev.off()  

dfcorrs<-enc_rfres[enc_rfres$condition%in%best_kds,]

#dfcorrs<-rfres
dfcorrs<-unique(dfcorrs[,c("RBP","cell","method","xy","feature","corrp_tst")])
dfcorrs<-dfcorrs[dfcorrs$feature=="RNA_log2FC",]
dfcorrs<-dfcorrs[dfcorrs$xy=="all",]

#dfcorrs<-dfcorrs[dfcorrs$method=="shRNA",]
#dfcorrs<-dfcorrs[dfcorrs$cell=="K562",]
dfcorrs<-dfcorrs[order(dfcorrs$corrp_tst,decreasing = T),]
dfcorrs$type_exp<-interaction(dfcorrs$cell,dfcorrs$method)
ago<-aggregate(dfcorrs$corrp_tst,by=list(dfcorrs$RBP),mean)
dfcorrs_cp<-dfcorrs
aribipisi<-ago[,1][order(ago[,2],decreasing = T)]
dfcorrs$RBP<-factor(dfcorrs$RBP,levels = aribipisi)
dfcorrs$method<-interaction(dfcorrs$cell,dfcorrs$method)
dfcorrs<-aggregate(dfcorrs$corrp_tst,by=list(dfcorrs$RBP,dfcorrs$method),mean)
colnames(dfcorrs)<-c("RBP","method","corrp_tst")
kdos2<-kdos
names(kdos2)<-sapply(strsplit(names(kdos2),"_version"),"[[",1)
dfcorrs$log2fc<-kdos2[paste(dfcorrs$RBP,sapply(strsplit(as.character(dfcorrs$method),"[.]"),FUN = function(x){paste(x[1],x[2],sep="_")}),sep = "_")]
#colsi<-setNames(object = alpha(c("black","blue","purple","red","orange"),.8),nm = c("all","x","-xy","y","xy"))
a<-ggplot(dfcorrs,aes(x=corrp_tst,fill=method)) + geom_density()
a<-a + theme_classic() +
  xlab("correlation with test data") +
  ylab("Density") +
  theme(axis.title.x = element_text(size=24),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_fill_manual(values = alpha(c("firebrick1","blue","dark red","dark blue"),.8),"method")
pdf(file = "Figures/Supplementary_Figure_6B.pdf",width = 7,height = 5)
#a+facet_wrap(feature~.,nrow = 2)
a
dev.off()

a<-ggplot(dfcorrs,aes(x=corrp_tst)) + geom_histogram(bins=40)
a<-a + theme_classic() +
  xlab("Pearson correlation with test data") +
  ylab("Datasets") +
  theme(axis.title.x = element_text(size=24),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  
#scale_fill_manual(values = alpha(c("firebrick1","blue","dark red","dark blue"),.8),"method")
pdf(file = "Figures/Fig4B.pdf",width = 7,height = 5)
#a+facet_wrap(feature~.,nrow = 2)
a
dev.off()


a<-ggplot(dfcorrs,aes(x=log2fc,fill=method)) + geom_density()
a<-a + theme_classic() +
  xlab("KD efficiency (log2FC)") +
  ylab("Density") +
  theme(axis.title.x = element_text(size=24),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_fill_manual(values = alpha(c("firebrick1","blue","dark red","dark blue"),.8),"method")


pdf(file = "Figures/Supplementary_Figure_6A.pdf",width = 7,height = 5)
#a+facet_wrap(feature~.,nrow = 2)
a
dev.off()


#check if control samples impact results

# 
# dfdGC<-dfd[dfd$Var1=="GCpct_cds",]
# dffGC<-dff[dfdGC$study]
# dfdGC$controls<-as.factor(as.character(paste(dffGC,collapse = ";")))
# #dfdGC<-dfdGC[dfdGC$RBP%in%names(table(dfdGC$RBP)[which(table(dfdGC$RBP)>1)]),]
# dfdGC<-dfdGC[order(dfdGC$importance,decreasing = T),]
# dfdGC$controls<-as.character(dfdGC$controls)
# dfdGC$controls<-factor(dfdGC$controls,levels=unique(dfdGC$controls))
# 
# boxplot(dfdGC$importance~dfdGC$controls)
# 
# anno<-anova(lm(importance ~ RBP + controls ,data = dfdGC),lm(importance ~ controls ,data = dfdGC),lm(importance ~ RBP ,data = dfdGC),test="Chisq")
# anno



#check metatx encore


#cell cycle
#newpap
tibook<-read.table("data/elife_supp.csv",header = F,sep = "\t",stringsAsFactors = F)
colnames(tibook)<-c("gene_name","m0","mu","gamma","ton","halflife","plateau")
tibook$gene_name<-sapply(strsplit(tibook$gene_name,"_"),"[[",1)
aa<-gcqnt_orig
tibo<-read.table("data/gene_classes.csv",header = T,skip = 3,sep = "\t",stringsAsFactors = F)
tibo<-tibo[,2:4]
colnames(tibo)<-c("gene_name","time_toneg","group")
tibo$gene_name<-sapply(strsplit(tibo$gene_name,"_"),"[[",1)
tibo<-tibo[tibo$gene_name%in%aa$gene_name,]
tibo$GC_cds<-aa$GCpct_cds[match(tibo$gene_name,aa$gene_name)]*100
tibook<-tibook[tibook$gene_name%in%aa$gene_name,]
tibook$GC_cds<-aa$GCpct_cds[match(tibook$gene_name,aa$gene_name)]*100

#gcqnt<-setNames(cut(tibo$GC_cds,breaks = quantile(tibo$GC_cds,prob=c(seq(0,1,length.out = 5))),include.lowest = T,right = T),nm = tibo$gene_name)
gcqnt<-setNames(cut(tibo$GC_cds,breaks = 3 ,include.lowest = T,right = T),nm = tibo$gene_name)
gcqnt<-gsub(gsub(gsub(gsub(as.character(gcqnt),pattern = "[(]",replacement = ""),pattern = "\\[",replacement = ""),pattern = "\\]",replacement = ""),pattern = ",",replacement = "-")
gcqnt<-factor(gcqnt)
gcqnt<-setNames(gcqnt,nm = tibo$gene_name)
tibo$gcqnt<-gcqnt[tibo$gene_name]
tibo$gcqnt<-as.character(tibo$gcqnt)
#boxplot(tibo$time_toneg~tibo$gcqnt)
#boxplot(tibo$time_toneg~cut(tibo$GC_cds,breaks = 4))
tibo2<-tibo
tibo2$gcqnt<-"all"
tibo<-rbind(tibo,tibo2)
tibo$gcqnt<-factor(tibo$gcqnt,levels=c("all",sort(unique(tibo$gcqnt))[-which(sort(unique(tibo$gcqnt))=="all")]))
collxy<-c("dark grey",colorRampPalette(colors = c("forestgreen","grey","purple"))(length(levels(tibo$gcqnt))-1))
names(collxy)<-levels(tibo$gcqnt)

btps<-unique(tibo$gcqnt)
comps<-unname(split(cbind(rep("all",length(btps)-1),as.character(btps[-4])),f = 1:(length(btps)-1)))
comps<-comps[c(2,1,3)]

imaging_decay<-ggplot(tibo,aes(y=time_toneg,x=gcqnt,color=gcqnt,group=gcqnt)) +stat_summary() + geom_boxplot()+
  geom_jitter(alpha=.3,size=.2) +
  #stat_summary(position=position_dodge(.9),size=.7) +
  ylab("Time to max mRNA decay") +
  xlab("%GC in CDS") +
  scale_color_manual(values = collxy,"%GC in CDS") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
#stat_compare_means(label="p.signif",method = "wilcox.test",comparisons = comps,label.y = 2.5 + seq(0,1,length.out = 5))
#facet_wrap(.~method) + coord_cartesian(ylim = c(0,1))
imaging_decay<-imaging_decay + stat_compare_means(label="p.signif",method = "wilcox.test",comparisons = comps)
#pdf(file = "Figures/Fig3_decayimaging.pdf",width = 9,height = 5)
#imaging_decay
#dev.off()




tibook<-read.table("data/elife_supp.csv",header = F,sep = "\t",stringsAsFactors = F)
colnames(tibook)<-c("gene_name","m0","mu","gamma","ton","halflife","plateau")
tibook$gene_name<-sapply(strsplit(tibook$gene_name,"_"),"[[",1)
aa<-gcqnt_orig
tibo<-read.table("data/gene_classes.csv",header = T,skip = 3,sep = "\t",stringsAsFactors = F)
tibo<-tibo[,2:4]
colnames(tibo)<-c("gene_name","time_toneg","group")
tibo$gene_name<-sapply(strsplit(tibo$gene_name,"_"),"[[",1)
tibo<-tibo[tibo$gene_name%in%aa$gene_name,]
tibo$GC_cds<-aa$GCpct_cds[match(tibo$gene_name,aa$gene_name)]*100
tibook<-tibook[tibook$gene_name%in%aa$gene_name,]
tibook$GC_cds<-aa$GCpct_cds[match(tibook$gene_name,aa$gene_name)]*100
#gcqnt<-setNames(cut(tibook$GC_cds,breaks = quantile(tibook$GC_cds,prob=c(seq(0,1,length.out = 5))),include.lowest = T,right = T),nm = tibook$gene_name)
gcqnt<-setNames(cut(tibook$GC_cds,breaks = 3 ,include.lowest = T,right = T),nm = tibook$gene_name)
gcqnt<-gsub(gsub(gsub(gsub(as.character(gcqnt),pattern = "[(]",replacement = ""),pattern = "\\[",replacement = ""),pattern = "\\]",replacement = ""),pattern = ",",replacement = "-")
gcqnt<-factor(gcqnt)
gcqnt<-setNames(gcqnt,nm = tibook$gene_name)


tibo$gcqnt<-gcqnt[tibo$gene_name]
tibo$gcqnt<-as.character(tibo$gcqnt)
#boxplot(tibo$time_toneg~tibo$gcqnt)
#boxplot(tibo$time_toneg~cut(tibo$GC_cds,breaks = 4))
#tibo<-rbind(tibo,tibo2)
tibo$gcqnt<-factor(tibo$gcqnt,levels=sort(unique(tibo$gcqnt)))
collxy<-colorRampPalette(colors = c("forestgreen","grey","purple"))(length(levels(tibo$gcqnt)))
names(collxy)<-levels(tibo$gcqnt)

btps<-unique(tibo$gcqnt)
comps<-list()
comps[[1]]<-c("47.9-61.1","34.6-47.9")
comps[[2]]<-c("61.1-74.4","47.9-61.1")


imaging_decay<-ggplot(tibo,aes(y=time_toneg,x=gcqnt,color=gcqnt,group=gcqnt)) +stat_summary() + geom_boxplot()+
  geom_jitter(alpha=.5,size=.2) +
  #stat_summary(position=position_dodge(.9),size=.7) +
  ylab("Time to min mRNA levels") +
  xlab("%GC in CDS") +
  scale_color_manual(values = collxy,"%GC in CDS") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
#stat_compare_means(label="p.signif",method = "wilcox.test",comparisons = comps,label.y = 2.5 + seq(0,1,length.out = 5))
#facet_wrap(.~method) + coord_cartesian(ylim = c(0,1))
imaging_decay<-imaging_decay + stat_compare_means(method = "wilcox.test",comparisons = comps,method.args = list(alternative="l"))
#pdf(file = "Figures/Fig3_decayimaging.pdf",width = 9,height = 5)
imaging_decaymin<-imaging_decay +  theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank())
#dev.off()


tibook$gcqnt<-gcqnt[tibook$gene_name]
tibook$gcqnt<-as.character(tibook$gcqnt)
#boxplot(tibook$time_toneg~tibook$gcqnt)
#boxplot(tibook$time_toneg~cut(tibook$GC_cds,breaks = 4))
tibook$gcqnt<-factor(tibook$gcqnt,levels=sort(unique(tibook$gcqnt)))

imaging_decay<-ggplot(tibook,aes(y=gamma,x=gcqnt,color=gcqnt,group=gcqnt))  + geom_boxplot() +
  geom_jitter(alpha=.5,size=.2) +
  #stat_summary(position=position_dodge(.9),size=.7) +
  ylab("Decay rate") +
  xlab("%GC in CDS") +
  scale_color_manual(values = collxy,"%GC in CDS") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
#stat_compare_means(method = "wilcox.test",comparisons = comps,label.y = 2.5 + seq(0,1,length.out = 5))
#facet_wrap(.~method) + coord_cartesian(ylim = c(0,1))
imaging_decay<-imaging_decay 
#pdf(file = "Figures/Fig3_decayimaging.pdf",width = 9,height = 5)
imaging_decaygamma<-imaging_decay + coord_cartesian(ylim=c(0,.15))+ stat_compare_means(method = "wilcox.test",comparisons = comps,label.y = mean(tibook$gamma,na.rm=T)*2 + seq(mean(tibook$gamma,na.rm=T)/2,mean(tibook$gamma,na.rm=T)*2,length.out = 4),method.args = list(alternative="g")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#dev.off()


imaging_decayhalf<-ggplot(tibook,aes(y=halflife,x=gcqnt,color=gcqnt,group=gcqnt))  + geom_boxplot() +
  geom_jitter(alpha=.5,size=.2) +
  #stat_summary(position=position_dodge(.9),size=.7) +
  ylab("Half-life") +
  xlab("%GC in CDS") +
  scale_color_manual(values = collxy,"%GC in CDS") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
#stat_compare_means(method = "wilcox.test",comparisons = comps,label.y = 2.5 + seq(0,1,length.out = 5))
#facet_wrap(.~method) + coord_cartesian(ylim = c(0,1))
imaging_decayhalf<-imaging_decayhalf + stat_compare_means(method = "wilcox.test",comparisons = comps,label.y = mean(tibook$halflife,na.rm=T)*2 + seq(mean(tibook$halflife,na.rm=T)/2,mean(tibook$halflife,na.rm=T)*2,length.out = 4),method.args = list(alternative="l")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#pdf(file = "Figures/Fig3_decayimaging.pdf",width = 9,height = 5)

#dev.off()

imaging_decay<-ggplot(tibook,aes(y=ton,x=gcqnt,color=gcqnt,group=gcqnt))  + geom_boxplot() +
  geom_jitter(alpha=.5,size=.2) +
  #stat_summary(position=position_dodge(.9),size=.7) +
  ylab("Onset of decay") +
  xlab("%GC in CDS") +
  scale_color_manual(values = collxy,"%GC in CDS") +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
#stat_compare_means(method = "wilcox.test",comparisons = comps,label.y = 2.5 + seq(0,1,length.out = 5))
#facet_wrap(.~method) + coord_cartesian(ylim = c(0,1))
imaging_decay<-imaging_decay + stat_compare_means(method = "wilcox.test",comparisons = comps,label.y = mean(tibook$ton,na.rm=T)*2 + seq(mean(tibook$ton,na.rm=T)/2,mean(tibook$ton,na.rm=T)*2,length.out = 4),method.args = list(alternative="l"))
#pdf(file = "Figures/Fig3_decayimaging.pdf",width = 9,height = 5)
imaging_decayton<- imaging_decay +  theme(axis.title.x=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank())
#dev.off()
pdf(file = "Figures/Fig4E.pdf",width = 6,height = 15)
plot_grid(imaging_decaymin,imaging_decaygamma,imaging_decayhalf,imaging_decayton,align = "v",ncol = 1)
dev.off()


#mouse

summary_multiDE_mouse<-get(load("data/mouseAllnoHet_multiDE_results_summary.RData"))
summary_multiDE_mouse$uniq$het<-NULL
plotti<-scatterplot_RiboRNAIntronExon(summary_multiDE = summary_multiDE_mouse)
dfa<-plotti$RiboRNA_tx_type$data
dfa$gene_type<-ifelse(test = dfa$single_exon,yes = "single_exon",no = "multi_exon")
colsi<-alpha(c("blue","dark red","gray43","gray10","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_NS")
dfa$GCpct_cds<-dfa$GCpct_cds*100
dfa$sig=0
dfa$sig[which(dfa$RNA_padj<.01 | dfa$RiboRNA_padj<.05)]=1
a<-ggplot(dfa,aes(x=RNA_log2FC,y=Ribo_log2FC,color=GCpct_cds,size=RiboRNA_padj,alpha=sig)) + geom_point()
a<-a + theme_bw() +
  ylab(paste("Ribo-seq "," log2FC",sep = "")) +
  xlab(paste("RNA-seq log2FC",sep = "")) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_gradient2(low = "forestgreen",mid = "grey",high = "purple",midpoint = mean(dfa$GCpct_cds),"%GC in CDS")+
  scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv,
                                                                               domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))
axm=2.5
a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + coord_cartesian(xlim = c(-axm,axm),ylim=c(-axm,axm))
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1)
pdf(file = "Figures/Fig6A.pdf",width = 9,height = 5)
ariborna + scale_alpha(guide="none")
dev.off()


a<-ggplot(dfa,aes(x=RNA_log2FC,y=Ribo_log2FC,color=tx_type_RiboRNA,size=RiboRNA_padj,shape=gene_type)) + geom_point()
a<-a + theme_bw() +
  ylab(paste("Ribo-seq "," log2FC",sep = "")) +
  xlab(paste("RNA-seq log2FC",sep = "")) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = colsi,"Tx_class")+
  scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv,
                                                                               domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))
axm=2.5
a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + coord_cartesian(xlim = c(-axm,axm),ylim=c(-axm,axm))
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1)
#pdf(file = "Figures/Fig4_RiboRNA_type.pdf",width = 9,height = 5)
#ariborna
#dev.off()

summary_multiDE_mouse$uniq$het<-NULL
plotti<-maplot_multiDE(summary_multiDE = summary_multiDE_mouse,assayType = "RNA",positional = T)

#pdf(file = "Figures/Supp_Fig4_MAposRNA.pdf",width = 9,height = 5)
#plotti
#dev.off()

#plotti<-maplot_multiDE(summary_multiDE = summary_multiDE_mouse,assayType = "Ribo",positional = T)
#pdf(file = "Figures/Supp_Fig4_MAposRibo.pdf",width = 9,height = 5)
#plotti
#dev.off()

load("data/rleall_stps.RData")
# add step to show violin plots (and add numbers) using .15 .1 and .2

ddx<-summary_multiDE$uniq$`48h`
ddx$tx_type_RiboRNA<-"mixed"
ddx$tx_type_RiboRNA[ddx$RNA_padj<.01 & ddx$RNA_log2FC<0 & ddx$Ribo_log2FC<0]<-"Concordant_down"
ddx$tx_type_RiboRNA[ddx$RNA_padj<.01 & ddx$RNA_log2FC>0 & ddx$Ribo_log2FC>0]<-"Concordant_up"
ddx$tx_type_RiboRNA[ddx$RiboRNA_padj<.05 & ddx$RiboRNA_log2FC<0]<-"TE_down"
ddx$tx_type_RiboRNA[ddx$RiboRNA_padj<.05 & ddx$RiboRNA_log2FC>0]<-"TE_up"
ddx$tx_type_RiboRNA[which((ddx$Ribo_padj>.2 | is.na(ddx$RNA_padj)) & (ddx$RNA_padj>.2| is.na(ddx$Ribo_padj))) ]<-"unchanging"
ddx$totake<-F
set.seed(667)
totak<-sample(which(ddx$tx_type_RiboRNA=="unchanging" & ddx$TPM_RNA>10 & ddx$TPM_Ribo>10 & !is.na(ddx$RNA_padj) & !is.na(ddx$Ribo_padj)),size = 1500,replace = F)
ddx$totake[totak]<-T

ordd<-ddx[order(ddx$RNA_padj),]
gnsup<-ordd[which(ordd$RNA_log2FC>0 & !is.na(ordd$RNA_padj) & !is.na(ordd$Ribo_padj)),"gene_name"][1:250]
gnsdown<-ordd[which(ordd$RNA_log2FC<0 & !is.na(ordd$RNA_padj) & !is.na(ordd$Ribo_padj)) ,"gene_name"][1:250]
ddx$type_toRle<-NA
ddx$type_toRle[ddx$totake]<-"unchanging"
ddx$type_toRle[ddx$gene_name%in%gnsup]<-"stabilized"
ddx$type_toRle[ddx$gene_name%in%gnsdown]<-"degraded"

ddx$totake[ddx$gene_name%in%c(gnsup,gnsdown)]<-T
ddx_hum<-ddx

df<-do.call(lapply(runlsok,function(x){
  
  mino<-min(x$startp_rna,na.rm = T)
  topad<-(mino+500+249)-dim(x$profilesavg)[1]
  
  if(topad<=0){
    rnapro<-x$profilesavg[(mino+500-250):(mino+500+249),1:6]
    ribopro<-x$profilesavg[(mino+500-250):(mino+500+249),7:12]
  }
  if(topad>0){
    rnapro<-rbind(x$profilesavg[(mino+500-250):dim(x$profilesavg)[1],1:6],matrix(NA,nrow = topad,ncol = 6))
    ribopro<-rbind(x$profilesavg[(mino+500-250):dim(x$profilesavg)[1],7:12],matrix(NA,nrow = topad,ncol = 6))
    
  }
  
  colnames(rnapro)<-c("DMSO","4h","8h","12h","24h","48h")
  colnames(ribopro)<-c("DMSO","4h","8h","12h","24h","48h")
  
  dfrn<-suppressMessages(melt(rnapro))
  dfrn$type<-"RNA"
  dfrn2<-suppressMessages(melt(ribopro))
  dfrn2$type<-"Ribo"
  df<-rbind(dfrn,dfrn2)
  df$position<-rep(1:500,dim(df)[1]/500)
  df$region<-"start"
  
  dfst<-df
  
  rnapro<-x$profilesendavg[2:501,1:6]
  ribopro<-x$profilesendavg[2:501,7:12]
  
  colnames(rnapro)<-c("DMSO","4h","8h","12h","24h","48h")
  colnames(ribopro)<-c("DMSO","4h","8h","12h","24h","48h")
  
  dfrn<-suppressMessages(melt(rnapro))
  dfrn$type<-"RNA"
  dfrn2<-suppressMessages(melt(ribopro))
  dfrn2$type<-"Ribo"
  df<-rbind(dfrn,dfrn2)
  df$position<-rep(1:500,dim(df)[1]/500)
  df$region<-"end"
  df<-rbind(dfst,df)
  df
}),what = rbind)
colnames(df)<-c("position2","variable","value","type","position","region")
df$gene_name<-sapply(strsplit(rownames(df),"[.]"),"[[",1)
df$tx_class<-setNames(ddx$type_toRle,ddx$gene_name)[df$gene_name]

#names(head(sort(anovRNA,decreasing = F),100))[ddx$type_toRle[match(names(head(sort(anovRNA,decreasing = F),100)),ddx$gene_name)]=="stabilized"]

diststart<-sapply(runlsok,function(x){median(x$pospointone[,1])-min(x$pospointone[,1])})
colss<-alpha(colorRampPalette(colors = c("dark grey","steelblue","dark blue"))(6),.8)
dfok<-df

agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
agg2<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),sd,na.rm=T)
dfok<-agg1
dfok$sd<-agg2$x
colnames(dfok)<-c("position","variable","tx_class","type","region","value","sd")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
pl<-ggplot(dfok,aes(x=position,y=value,color=variable)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  geom_hline(yintercept = 0,linetype=2)+
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

#pdf(file = "Figures/Fig5_aggrehuman.pdf",width = 14,height = 5)
#pl + facet_wrap(region~tx_class)  + scale_color_manual(values = colss,"Time after\ndegron activation") 
#dev.off()


dfok<-df
#dfok$position[dfok$region=="start"]<-dfok$position[dfok$region=="start"]-251
#dfok$position[dfok$region=="end"]<-dfok$position[dfok$region=="end"]-252
dfok<-dfok[dfok$gene_name=="CSRNP2",]
stpnts<-runlsok$CSRNP2$pospointone[,1]-min(runlsok$CSRNP2$startp_rna,na.rm = T)+250

datap<-data.frame(position=stpnts,value=0,variable=levels(dfok$variable))
agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
agg2<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),sd,na.rm=T)
dfok<-agg1
dfok<-agg1
dfok$sd<-agg2$x
colnames(dfok)<-c("position","variable","tx_class","type","region","value","sd")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
pl<-ggplot(dfok,aes(x=position,y=value,color=variable)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

pdf(file = "Figures/Supplementary_Figure_8.pdf",width = 14,height = 5)
pl + facet_wrap(region~.)  + scale_color_manual(values = colss,"Time after\ndegron activation")
dev.off()


dfok2<-dfok[dfok$region=="start",]
stpnts<-runlsok$CSRNP2$pospointone[,1]-min(runlsok$CSRNP2$startp_rna,na.rm = T)+250
#datap<-data.frame(position=stpnts,value=0,variable=levels(dfok$variable))
#datap<-rbind(datap,data.frame(position=stpnts,value=.15,variable=levels(dfok$variable)))
datap<-data.frame(position=stpnts,value=.15,variable=levels(dfok$variable))
datap$reach=.9
dfok2$reach<-.9
for( i in levels(dfok2$variable)){
  valsi<-dfok2$value[dfok2$variable==i]
  rici<-dfok2$reach[dfok2$variable==i]
  rici[which(valsi>.15)[1]:length(valsi)]<-.2
  dfok2$reach[dfok2$variable==i]<-rici
}

pl<-ggplot(dfok2,aes(x=position,y=value,color=variable,alpha=reach)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  scale_alpha_continuous(range=c(.2,.9))

pdf(file = "Figures/Fig5A3A.pdf",width = 8,height = 5)
pl + scale_color_manual(values = colss,"Time after\ndegron activation")+
  geom_vline(xintercept = stpnts,color=c(colss),linetype=2) + geom_point(data = datap,size=5,color=c(colss)) +
  scale_x_continuous(breaks = c(100,250,400),labels = c(-150,0,150))
dev.off()


reslm<-lapply(runlsok,function(x){
  x<-data.frame(x$pospointone)
  x$time<-c(0,4,8,16,24,48)
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time,data = x)$coefficients[2])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time,data = x)$coefficients[2])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  rescoef
  
})

coefsRNA<-sapply(reslm,"[[",1)
coefsRibo<-sapply(reslm,"[[",2)
anovRNA<-sapply(reslm,"[[",3)
anovRibo<-sapply(reslm,"[[",4)
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
#boxplot(coefsRNA~ginsclass,ylim=c(-20,20))
#wilcox.test(coefsRNA[ginsclass=="degraded"],coefsRNA[ginsclass=="unchanging"])$p.value
#wilcox.test(coefsRNA[ginsclass=="stabilized"],coefsRNA[ginsclass=="unchanging"])$p.value
stpoint<-sapply(runlsok,function(x){x$startp_rna})
stpoint<-sort(stpoint,decreasing = T)
ginsclass<-ddx$type_toRle[match(names(stpoint),ddx$gene_name)]
gcclass<-ddx$GCpct_cds[match(names(stpoint),ddx$gene_name)]
fcclass<-ddx$RNA_log2FC[match(names(stpoint),ddx$gene_name)]

#boxplot(stpoint~ginsclass,ylim=c(0,600))
#wilcox.test(stpoint[ginsclass=="degraded"],stpoint[ginsclass=="unchanging"])$p.value
#wilcox.test(stpoint[ginsclass=="stabilized"],stpoint[ginsclass=="unchanging"])$p.value

x<-runlsok$CSRNP2
x<-data.frame(x$pospointone)
x$time<-c(0,4,8,16,24,48)
colnames(x)<-c("RNA","Ribo","time")
#x$RNA<-as.numeric(scale(x$RNA))
pl<-ggplot(x,aes(x=time,y=RNA)) + geom_point(size=3) + geom_smooth(method='lm',se = F) +
  xlab("hours after degron induction") +
  ylab("coverage starting position") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_regline_equation(label.y = min(x$RNA), aes(label = ..eq.label..),size=9) 

pdf(file = "Figures/Fig5A4A_lm.pdf",width = 10,height = 5)
pl + scale_x_continuous(breaks=c(0,4,8,16,24,48))
dev.off()


coefsRNA<-sapply(reslm,"[[",1)
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("degraded", "unchanging") )
colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(coverage start position)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = -1.5+1:3,method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


pdf(file = "Figures/Fig5A5A.pdf",width = 10,height = 5)
pl + coord_cartesian(ylim = c(-6,4))
dev.off()

stpoint<-sapply(runlsok,function(x){x$startp_rna})
stpoint<-sort(stpoint,decreasing = T)
ginsclass<-ddx$type_toRle[match(names(stpoint),ddx$gene_name)]

dfok<-cbind.data.frame(stpoint,ginsclass)
colnames(dfok)<-c("star_coverage","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
#my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("unchanging", "degraded") )

colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,star_coverage,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("start coverage\nposition (nt)") +
  scale_fill_manual(values = colsi2) +
  #geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 500+c(100,150,200),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))


#pdf(file = "Figures/Fig5A5B_startcov.pdf",width = 10,height = 5)
#pl + coord_cartesian(ylim = c(0,900))
#dev.off()

# in addition

meancovst<-lapply(runlsok, function(x){
  topad<-300-(dim(x$profilesavg)[1]-(x$startp_rna+500))
  if(topad>0){x$profilesavg<-rbind(x$profilesavg,matrix(0,nrow = topad,ncol=ncol(x$profilesavg)))}
  #x$profilesavg<-apply(x$profilesavg,2,function(x){x/max(x,na.rm = T)})
  y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
  
  y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
  x<-data.frame(y)
  x$time<-c(0,4,8,16,24,48)
  
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(x[1,1]==0){x[1,1]<-x[2,1]}
  if(x[1,2]==0){x[1,2]<-x[2,2]}
  x[,1]<-log2(x[,1]/x[1,1])
  x[,2]<-log2(x[,2]/x[1,2])
  if(length(which(is.infinite(x[,1])))>0){x[which(is.infinite(x[,1])),1]<-0}
  if(length(which(is.infinite(x[,2])))>0){x[which(is.infinite(x[,2])),2]<-0}
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time -1,data = x)$coefficients[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time -1,data = x)$coefficients[1])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  resl<-list(y,rescoef)
  names(resl)<-c("avg_cov","lmres")
  resl
  
})



coefsRNA<-sapply(meancovst,function(x){x$lmres[1]})
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
#my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("unchanging", "degraded") )
colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(log2FC coverage)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 0.02+c(.02,.03,.04),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))


pdf(file = "Figures/Fig5A5B.pdf",width = 10,height = 5)
#pl + coord_cartesian(ylim = c(-6,4))
pl
dev.off()



dfok<-df
#dfok$position[dfok$region=="start"]<-dfok$position[dfok$region=="start"]-251
#dfok$position[dfok$region=="end"]<-dfok$position[dfok$region=="end"]-252
dfok<-dfok[dfok$gene_name=="CSRNP2",]
stpnts<-runlsok$CSRNP2$pospointone[,1]-min(runlsok$CSRNP2$startp_rna,na.rm = T)+250

datap<-data.frame(position=stpnts,value=.15,variable=levels(dfok$variable))
agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
agg2<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),sd,na.rm=T)
dfok<-agg1
dfok<-agg1
dfok$sd<-agg2$x
colnames(dfok)<-c("position","variable","tx_class","type","region","value","sd")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
dfok2<-dfok[dfok$region=="start",]
pl<-ggplot(dfok2,aes(x=position,y=value,color=variable)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

#pdf(file = "Figures/Fig5_metaexample.pdf",width = 9,height = 5)
#pl + scale_color_manual(values = colss,"Time after\ndegron activation")+
  #geom_vline(xintercept = stpnts,color=c(colss),linetype=2) + geom_point(data = datap,size=5,color=c(colss))
#dev.off()



x<-runlsok$CSRNP2
y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
newpcov<-y[,1]
x<-data.frame(y)
x$time<-c(0,4,8,16,24,48)
colnames(x)<-c("RNA","Ribo","time")
x[,1]<-log2(x[,1]/x[1,1])
x[,2]<-log2(x[,2]/x[1,2])


#x$RNA<-as.numeric(scale(x$RNA))
pl<-ggplot(x,aes(x=time,y=RNA)) + geom_point(size=3) + geom_smooth(method='lm',formula = y ~ x -1,se = F) +
  xlab("hours after degron induction") +
  ylab("log2FC coverage\nstarting position") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_regline_equation(label.y = min(x$RNA), formula = y ~ x -1 , aes(label = ..eq.label..),size=9) 

pdf(file = "Figures/Fig5A4B_lm2.pdf",width = 10,height = 5)
pl + scale_x_continuous(breaks=c(0,4,8,16,24,48))
dev.off()




dfok<-df
#dfok$position[dfok$region=="start"]<-dfok$position[dfok$region=="start"]-251
#dfok$position[dfok$region=="end"]<-dfok$position[dfok$region=="end"]-252
dfok<-dfok[dfok$gene_name=="CSRNP2",]

datap<-data.frame(position=0,value=newpcov,variable=levels(dfok$variable))
agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
agg2<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),sd,na.rm=T)
dfok<-agg1
dfok<-agg1
dfok$sd<-agg2$x
colnames(dfok)<-c("position","variable","tx_class","type","region","value","sd")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
dfok2<-dfok[dfok$region=="start",]

datap$reach=.9
dfok2$reach<-.3
for( i in levels(dfok2$variable)){
  rici<-dfok2$reach[dfok2$variable==i]
  rici[100:400]<-.9
  dfok2$reach[dfok2$variable==i]<-rici
}

pl<-ggplot(dfok2,aes(x=position,y=value,color=variable,alpha=reach)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position centered on coverage start") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  geom_hline(yintercept = newpcov,color=c(colss),linetype=2)+
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

pdf(file = "Figures/Fig5A3B.pdf",width = 8,height = 5)
pl + scale_color_manual(values = colss,"Time after\ndegron activation")+
  geom_point(data = datap,size=5,color=c(colss)) + scale_x_continuous(breaks = c(100,250,400),labels = c(-150,0,150)) +
  geom_vline(xintercept = c(100,400))
dev.off()

pl<-ggplot(dfok2,aes(x=position,y=value,color=variable)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position centered on coverage start") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  #geom_hline(yintercept = newpcov,color=c(colss),linetype=2)+
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

pdf(file = "Figures/Fig5A2.pdf",width = 8,height = 5)
pl + scale_color_manual(values = colss,"Time after\ndegron activation")+
  scale_x_continuous(breaks = c(100,250,400),labels = c(-150,0,150)) 

dev.off()




dfok<-df
#dfok$position[dfok$region=="start"]<-dfok$position[dfok$region=="start"]-251
#dfok$position[dfok$region=="end"]<-dfok$position[dfok$region=="end"]-252
dfok<-dfok[dfok$gene_name=="CSRNP2",]

datap<-data.frame(position=0,value=newpcov,variable=levels(dfok$variable))
agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
dfok<-agg1
colnames(dfok)<-c("position","tx_class","type","region","value")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
dfok2<-dfok[dfok$region=="start",]


pl<-ggplot(dfok2,aes(x=position,y=value,color="black")) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position centered on coverage start") +
  scale_color_manual(values="black","Time after\ndegron activation") +
  theme_classic() +
  #geom_hline(yintercept = newpcov,color=c(colss),linetype=2)+
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

pdf(file = "Figures/Fig5A1.pdf",width = 8,height = 5)
pl + scale_x_continuous(breaks = c(100,250,400),labels = c(-150,0,150)) 

dev.off()



#ADD metatx mouse and encode

load("data/silver_rleall_stps.RData")


ddx<-get(load("data/mouseAllnoHet_multiDE_results_summary.RData"))
load("data/gencode.vM22.annotation.gtf_Rannot")
ddx<-ddx$uniq$cKO

ddx$tx_type_RiboRNA<-"mixed"
ddx$tx_type_RiboRNA[ddx$RNA_padj<.01 & ddx$RNA_log2FC<0 & ddx$Ribo_log2FC<0]<-"Concordant_down"
ddx$tx_type_RiboRNA[ddx$RNA_padj<.01 & ddx$RNA_log2FC>0 & ddx$Ribo_log2FC>0]<-"Concordant_up"
ddx$tx_type_RiboRNA[ddx$RiboRNA_padj<.05 & ddx$RiboRNA_log2FC<0]<-"TE_down"
ddx$tx_type_RiboRNA[ddx$RiboRNA_padj<.05 & ddx$RiboRNA_log2FC>0]<-"TE_up"
ddx$tx_type_RiboRNA[which((ddx$Ribo_padj>.2 | is.na(ddx$RNA_padj)) & (ddx$RNA_padj>.2| is.na(ddx$Ribo_padj))) ]<-"unchanging"


ddx$totake<-F
set.seed(667)
totak<-sample(which(ddx$tx_type_RiboRNA=="unchanging" & ddx$TPM_RNA>10 & ddx$TPM_Ribo>10 & !is.na(ddx$RNA_padj) & !is.na(ddx$Ribo_padj)),size = 1500,replace = F)
ddx$totake[totak]<-T

ordd<-ddx[order(ddx$RNA_padj),]
gnsup<-ordd[which(ordd$RNA_log2FC>0 & !is.na(ordd$RNA_padj) & !is.na(ordd$Ribo_padj)),"gene_name"][1:250]
gnsdown<-ordd[which(ordd$RNA_log2FC<0 & !is.na(ordd$RNA_padj) & !is.na(ordd$Ribo_padj)) ,"gene_name"][1:250]
ddx$type_toRle<-NA
ddx$type_toRle[ddx$totake]<-"unchanging"
ddx$type_toRle[ddx$gene_name%in%gnsup]<-"stabilized"
ddx$type_toRle[ddx$gene_name%in%gnsdown]<-"degraded"

ddx$totake[ddx$gene_name%in%c(gnsup,gnsdown)]<-T
df<-do.call(lapply(runlsok,function(x){
  
  mino<-min(x$startp_rna,na.rm = T)
  topad<-(mino+500+249)-dim(x$profilesavg)[1]
  cnms<-do.call(strsplit(colnames(x$profilesavg),";"),what = rbind)
  if(topad<=0){
    rnapro<-x$profilesavg[(mino+500-250):(mino+500+249),cnms[,1]=="RNA"]
    ribopro<-x$profilesavg[(mino+500-250):(mino+500+249),cnms[,1]=="Ribo"]
  }
  if(topad>0){
    rnapro<-rbind(x$profilesavg[(mino+500-250):dim(x$profilesavg)[1],cnms[,1]=="RNA"],matrix(NA,nrow = topad,ncol = sum(cnms[,1]=="RNA")))
    ribopro<-rbind(x$profilesavg[(mino+500-250):dim(x$profilesavg)[1],cnms[,1]=="Ribo"],matrix(NA,nrow = topad,ncol = sum(cnms[,1]=="Ribo")))
    
  }
  
  colnames(rnapro)<-cnms[cnms[,1]=="RNA",2]
  colnames(ribopro)<-cnms[cnms[,1]=="Ribo",2]
  
  dfrn<-suppressMessages(melt(rnapro))
  dfrn$type<-"RNA"
  dfrn2<-suppressMessages(melt(ribopro))
  dfrn2$type<-"Ribo"
  df<-rbind(dfrn,dfrn2)
  df$position<-rep(1:500,dim(df)[1]/500)
  df$region<-"start"
  
  dfst<-df
  
  rnapro<-x$profilesendavg[2:501,cnms[,1]=="RNA"]
  ribopro<-x$profilesendavg[2:501,cnms[,1]=="Ribo"]
  
  colnames(rnapro)<-cnms[cnms[,1]=="RNA",2]
  colnames(ribopro)<-cnms[cnms[,1]=="Ribo",2]
  
  dfrn<-suppressMessages(melt(rnapro))
  dfrn$type<-"RNA"
  dfrn2<-suppressMessages(melt(ribopro))
  dfrn2$type<-"Ribo"
  df<-rbind(dfrn,dfrn2)
  df$position<-rep(1:500,dim(df)[1]/500)
  df$region<-"end"
  df<-rbind(dfst,df)
  df
}),what = rbind)

colnames(df)<-c("position2","variable","value","type","position","region")
df$gene_name<-sapply(strsplit(rownames(df),"[.]"),"[[",1)
df$tx_class<-setNames(ddx$type_toRle,ddx$gene_name)[df$gene_name]

#names(head(sort(anovRNA,decreasing = F),100))[ddx$type_toRle[match(names(head(sort(anovRNA,decreasing = F),100)),ddx$gene_name)]=="stabilized"]

colss<-alpha(colorRampPalette(colors = c("dark grey","steelblue","dark blue"))(length(names(table(df$variable)))),.8)
dfok<-df

agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
agg2<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),sd,na.rm=T)
dfok<-agg1
dfok$sd<-agg2$x
colnames(dfok)<-c("position","variable","tx_class","type","region","value","sd")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
pl<-ggplot(dfok,aes(x=position,y=value,color=variable)) + geom_line(size=1.5) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  geom_hline(yintercept = 0,linetype=2)+
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
#pdf(file = "Figures/Fig5_aggremouse.pdf",width = 14,height = 5)
#pl + facet_wrap(region~tx_class)  + scale_color_manual(values = colss,"condition") 
#dev.off()


dfok<-df
#dfok$position[dfok$region=="start"]<-dfok$position[dfok$region=="start"]-251
#dfok$position[dfok$region=="end"]<-dfok$position[dfok$region=="end"]-252
dfok<-dfok[dfok$gene_name=="Ctxn1",]
stpnts<-runlsok$Ctxn1$pospointone[,1]-min(runlsok$Ctxn1$startp_rna,na.rm = T)+250

datap<-data.frame(position=stpnts,value=.15,variable=levels(dfok$variable))


x<-runlsok$Ctxn1
y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
newpcov<-y[,1]
datap2<-data.frame(position=0,value=newpcov,variable=levels(dfok$variable))


agg1<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),mean,na.rm=T)
agg2<-aggregate(dfok$value,by=list(dfok$position,dfok$variable,dfok$tx_class,dfok$type,dfok$region),sd,na.rm=T)
dfok<-agg1
dfok<-agg1
dfok$sd<-agg2$x
colnames(dfok)<-c("position","variable","tx_class","type","region","value","sd")
dfok<-dfok[dfok$type=="RNA",]
dfok$region<-factor(dfok$region,levels = c("start","end"))
pl<-ggplot(dfok,aes(x=position,y=value,color=variable)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

pdf(file = "Figures/Supplementary_Figure_11.pdf",width = 14,height = 5)
pl + facet_wrap(region~.)  + scale_color_manual(values = colss,"")
dev.off()


dfok2<-dfok[dfok$region=="start",]
pl<-ggplot(dfok2,aes(x=position,y=value,color=variable)) + geom_line(size=1) +
  ylab("normalized RNA-seq coverage") +
  xlab("position (nt)") +
  #scale_color_manual(values = colsi,"") +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))

pdf(file = "Figures/Fig6B.pdf",width = 8,height = 5)
pl + scale_color_manual(values = colss,"condition")+
  geom_vline(xintercept = stpnts,color=c(colss),linetype=2) + geom_point(data = datap,size=5,color=c(colss)) +
  geom_hline(yintercept = newpcov,color=c(colss),linetype=2) + geom_point(data = datap2,size=5,color=c(colss))
dev.off()
# get diststart thingy and try it

stpoint<-sapply(runlsok,function(x){x$startp_rna})
stpoint<-sort(stpoint,decreasing = T)
ginsclass<-ddx$type_toRle[match(names(stpoint),ddx$gene_name)]
#boxplot(coefsRNA~ginsclass,ylim=c(-20,20))
#wilcox.test(coefsRNA[ginsclass=="degraded"],coefsRNA[ginsclass=="unchanging"])$p.value
#wilcox.test(coefsRNA[ginsclass=="stabilized"],coefsRNA[ginsclass=="unchanging"])$p.value


coefsRNA<-sapply(runlsok,function(x){diff(x$pospointone[,1])})
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("degraded", "unchanging") )
colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(draw_quantiles = .5)+
  xlab("") +
  ylab("delta start coverage") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = -200+c(20,50,80),method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


pdf(file = "Figures/Supplementary_Figure_10A.pdf",width = 10,height = 5)
pl + coord_cartesian(ylim = c(-200,50))
dev.off()

stpoint<-sapply(runlsok,function(x){x$startp_rna})
stpoint<-sort(stpoint,decreasing = T)
ginsclass<-ddx$type_toRle[match(names(stpoint),ddx$gene_name)]

dfok<-cbind.data.frame(stpoint,ginsclass)
colnames(dfok)<-c("star_coverage","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))

colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,star_coverage,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("start coverage\nposition (nt)") +
  scale_fill_manual(values = colsi2) +
  #geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 300+c(100,150,200),method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


#pdf(file = "Figures/Fig5_suppstartcovmouse.pdf",width = 10,height = 5)
#pl + coord_cartesian(ylim = c(0,900))
#dev.off()


meancovst<-lapply(runlsok, function(x){
  topad<-300-(dim(x$profilesavg)[1]-(x$startp_rna+500))
  if(topad>0){x$profilesavg<-rbind(x$profilesavg,matrix(0,nrow = topad,ncol=ncol(x$profilesavg)))}
  x$profilesavg<-apply(x$profilesavg,2,function(x){x/max(x,na.rm = T)})
  y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
  
  y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
  log(y[2,1]/y[1,1])
  
})



coefsRNA<-sapply(meancovst,function(x){x})
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
#my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("unchanging", "degraded") )
colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("log2FC 5' coverage") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 1+c(.5,1,1.5),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))


pdf(file = "Figures/Supplementary_Figure_10B.pdf",width = 10,height = 5)
pl + coord_cartesian(ylim = c(-1,4))
dev.off()


#Encode


annotation_file<-"data/gencode.v32.annotation.gtf.gz_Rannot"
load(annotation_file)


load("data/enc_startcov")

pvs<-lapply(enc_startcov,function(x){
  x$diffp<-x$pospointone2-x$pospointone1
  pvst<-wilcox.test(x$log_cov[x$tx_type=="stabilized"],x$log_cov[x$tx_type=="unchanging"],alternative = "g")$p.value
  pvdeg<-wilcox.test(x$log_cov[x$tx_type=="unchanging"],x$log_cov[x$tx_type=="degraded"],alternative = "g")$p.value
  pvst2<-wilcox.test(x$diffp[x$tx_type=="stabilized"],x$diffp[x$tx_type=="unchanging"],alternative = "l")$p.value
  pvdeg2<-wilcox.test(x$diffp[x$tx_type=="unchanging"],x$diffp[x$tx_type=="degraded"],alternative = "l")$p.value
  
  c(pvst,pvdeg,pvst2,pvdeg2)
})

pvs<-do.call(pvs,what = rbind.data.frame)
colnames(pvs)<-c("pval_cov_stab","pval_cov_deg","pval_pos_stab","pval_pos_deg")
pvs$sample<-names(enc_startcov)

df<-dfd[dfd$Var1%in%c("GCpct_cds"),]
df$study<-paste(df$RBP,df$cell,df$method,sep = "_")
pvs$study<-sapply(strsplit(pvs$sample,"_"),function(x){paste(x[-4],collapse="_")})
pvsok<-pvs[pvs$sample%in%best_kds,]
pvsok$importance_GCcds<-setNames(df$importance,df$study)[pvsok$study]

pvsok$importance_GCpct_cds<-cut(pvsok$importance_GCcds,breaks = quantile(pvsok$importance_GCcds,prob=c(seq(0,1,length.out = 5))),include.lowest = T,right = T)
#pvsok$importance_GCpct_cds<-cut(pvsok$importance_GCcds,breaks = 4,include.lowest = T,right = T)
pvsok$importance_GCpct_cds<-gsub(gsub(gsub(gsub(as.character(pvsok$importance_GCpct_cds),pattern = "[(]",replacement = ""),pattern = "\\[",replacement = ""),pattern = "\\]",replacement = ""),pattern = ",",replacement = "-")
tbb<-names(table(pvsok$importance_GCpct_cds))
#my_comparisons <- list( c("0.0233-0.0454", "0.00111-0.0233"),c("0.0454-0.0674", "0.00111-0.0233"),c("0.0674-0.0896","0.00111-0.0233") )
my_comparisons <- list( c(tbb[1], tbb[2]),c(tbb[1], tbb[3]),c(tbb[1],tbb[4]))
colsi2<-colorRampPalette(c("dark grey","dark blue"))(length(table(pvsok$importance_GCpct_cds)))
pl<-ggplot(pvsok,aes(importance_GCpct_cds,pval_cov_stab,fill=importance_GCpct_cds)) + geom_violin(scale = "area",draw_quantiles = .5)+
  ylab("5' coverage differences\nstabilized vs. unchanging (p-value)") +
  xlab("Importance\n%GC in CDS") +
  scale_fill_manual(values = colsi2,"Importance\n%GC in CDS") +
  #geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 8+c(0,1,2),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))

pdf(file = "Figures/Fig6C.pdf",width = 9,height = 7)

pl + scale_y_continuous(trans = scales::trans_new(name = "test",transform = trans,inverse = inv,
                                                   domain = c(1e-100, Inf)),breaks = c(10e-11,10e-8,10e-5,10e-2)) 
dev.off()
#boxplot(-log10(pvsok$pval_cov_stab)~cut(pvsok$importance_GCcds,breaks = 4))
#boxplot(-log10(pvsok$pval_cov_stab)~cut(pvsok$importance_GCcds,breaks = quantile(pvsok$importance_GCcds,prob=c(seq(0,1,length.out = 5))),include.lowest = T,right = T))

#make boxplot and prepare 1 stabilized, one degraded, and one unchanged mRNA to show 5' and 3' coverage of

load("data/enc_profiles.RData")
resprofs_orig<-resprofs

resprofs<-resprofs_orig


#ex_stab="KRT19"
ex_stab="GBP2"
#ex_unc="ZNF317"
ex_unc="DNAJC14"
#ex_unc="DPM1"

ex_deg="WARS"

#ZNF317

rbps<-sapply(strsplit(names(resprofs),"_"),"[[",1)
cell<-sapply(strsplit(names(resprofs),"_"),"[[",2)
method<-sapply(strsplit(names(resprofs),"_"),"[[",3)
cell<-factor(cell,levels = c("K562","HepG2"))
method<-factor(cell,levels = c("shRNA","CRISPR"))
ord<-order(rbps,cell,method)
resprofs<-resprofs[ord]
resprofs<-resprofs[!duplicated(rbps[ord])]

rbps<-sapply(strsplit(names(resprofs),"_"),"[[",1)
cell<-sapply(strsplit(names(resprofs),"_"),"[[",2)
method<-sapply(strsplit(names(resprofs),"_"),"[[",3)



ex_stabprof<-do.call(lapply(resprofs,function(x){x$profiles[[ex_stab]][400:650,]}),what = cbind)
ex_degprof<-do.call(lapply(resprofs,function(x){x$profiles[[ex_deg]][752:1252,]}),what = cbind)


prof_start<-do.call(lapply(resprofs,function(x){x$profiles[[ex_stab]][470:720,]}),what = cbind)
prof_end<-do.call(lapply(resprofs,function(x){x$profilesend[[ex_stab]][750:1000,]}),what = cbind)


ctrls<-prof_start[,colnames(prof_start)=="ctrl"]
kds<-prof_start[,colnames(prof_start)=="kd"]
dups<-duplicated(t(ctrls))
ctrls<-ctrls[,!dups]


colnames(ctrls)<-paste(colnames(ctrls),1:length(colnames(ctrls)),sep = "_")
colnames(kds)<-rbps

impgc<-setNames(df$importance,df$study)[sapply(strsplit(names(resprofs),"_"),function(x){paste(x[-4],collapse="_")})]
ordgc<-order(impgc,decreasing = F)

kdsok<-kds[,ordgc]
set.seed(665)
ctrlsok<-ctrls
#ctrlsok<-ctrls[,sample(1:length(ctrls[1,]),size = 10,replace = F)]

set.seed(666)
impgc2<-impgc
#sampei<-sort(sample(1:length(kdsok[1,]),size = 40,replace = F))
#impgc2<-impgc[ordgc][sampei]
impgc2["ctrl"]<-0
names(impgc2)<-sapply(strsplit(names(impgc2),"_"),"[[",1)
#kdsok<-kdsok[,sampei]

ddf<-melt(cbind(ctrlsok,kdsok))
ddf$type<-"RBP"
ddf$type[grep(ddf$Var2,pattern="ctrl")]<-"control"
ddf$type<-factor(ddf$type,levels = c("RBP","control"))
ddf$type2<-sapply(strsplit(as.character(ddf$Var2),"_"),"[[",1)
ddf$GCimp<-impgc2[as.character(ddf$type2)]
ddf$region<-"start"
ddfs<-ddf

ctrls<-prof_end[,colnames(prof_end)=="ctrl"]
kds<-prof_end[,colnames(prof_end)=="kd"]
dups<-duplicated(t(ctrls))
ctrls<-ctrls[,!dups]


colnames(ctrls)<-paste(colnames(ctrls),1:length(colnames(ctrls)),sep = "_")
colnames(kds)<-rbps

impgc<-setNames(df$importance,df$study)[sapply(strsplit(names(resprofs),"_"),function(x){paste(x[-4],collapse="_")})]
ordgc<-order(impgc,decreasing = F)

kdsok<-kds[,ordgc]
set.seed(665)
ctrlsok<-ctrls
#ctrlsok<-ctrls[,sample(1:length(ctrls[1,]),size = 10,replace = F)]

set.seed(666)
impgc2<-impgc
#sampei<-sort(sample(1:length(kdsok[1,]),size = 40,replace = F))
#impgc2<-impgc[ordgc][sampei]
impgc2["ctrl"]<-0
names(impgc2)<-sapply(strsplit(names(impgc2),"_"),"[[",1)
#kdsok<-kdsok[,sampei]

ddf<-melt(cbind(ctrlsok,kdsok))
ddf$type<-"RBP"
ddf$type[grep(ddf$Var2,pattern="ctrl")]<-"control"
ddf$type<-factor(ddf$type,levels = c("RBP","control"))
ddf$type2<-sapply(strsplit(as.character(ddf$Var2),"_"),"[[",1)
ddf$GCimp<-impgc2[as.character(ddf$type2)]
ddf$region<-"end"
ddf<-rbind(ddfs,ddf)
ddf$region<-factor(ddf$region,levels=c("start","end"))
ddf$typetx<-paste("stabilized",ex_stab,sep="_")

ddfst<-ddf


prof_start<-do.call(lapply(resprofs,function(x){x$profiles[[ex_deg]][500:750,]}),what = cbind)
prof_end<-do.call(lapply(resprofs,function(x){x$profilesend[[ex_deg]][750:1000,]}),what = cbind)


ctrls<-prof_start[,colnames(prof_start)=="ctrl"]
kds<-prof_start[,colnames(prof_start)=="kd"]
dups<-duplicated(t(ctrls))
ctrls<-ctrls[,!dups]


colnames(ctrls)<-paste(colnames(ctrls),1:length(colnames(ctrls)),sep = "_")
colnames(kds)<-rbps

impgc<-setNames(df$importance,df$study)[sapply(strsplit(names(resprofs),"_"),function(x){paste(x[-4],collapse="_")})]
ordgc<-order(impgc,decreasing = F)

kdsok<-kds[,ordgc]
set.seed(665)
ctrlsok<-ctrls
#ctrlsok<-ctrls[,sample(1:length(ctrls[1,]),size = 10,replace = F)]

set.seed(666)
impgc2<-impgc
#sampei<-sort(sample(1:length(kdsok[1,]),size = 40,replace = F))
#impgc2<-impgc[ordgc][sampei]
impgc2["ctrl"]<-0
names(impgc2)<-sapply(strsplit(names(impgc2),"_"),"[[",1)
#kdsok<-kdsok[,sampei]

ddf<-melt(cbind(ctrlsok,kdsok))
ddf$type<-"RBP"
ddf$type[grep(ddf$Var2,pattern="ctrl")]<-"control"
ddf$type<-factor(ddf$type,levels = c("RBP","control"))
ddf$type2<-sapply(strsplit(as.character(ddf$Var2),"_"),"[[",1)
ddf$GCimp<-impgc2[as.character(ddf$type2)]
ddf$region<-"start"
ddfs<-ddf

ctrls<-prof_end[,colnames(prof_end)=="ctrl"]
kds<-prof_end[,colnames(prof_end)=="kd"]
dups<-duplicated(t(ctrls))
ctrls<-ctrls[,!dups]


colnames(ctrls)<-paste(colnames(ctrls),1:length(colnames(ctrls)),sep = "_")
colnames(kds)<-rbps

impgc<-setNames(df$importance,df$study)[sapply(strsplit(names(resprofs),"_"),function(x){paste(x[-4],collapse="_")})]
ordgc<-order(impgc,decreasing = F)

kdsok<-kds[,ordgc]
set.seed(665)
ctrlsok<-ctrls
#ctrlsok<-ctrls[,sample(1:length(ctrls[1,]),size = 10,replace = F)]

set.seed(666)
impgc2<-impgc
#sampei<-sort(sample(1:length(kdsok[1,]),size = 40,replace = F))
#impgc2<-impgc[ordgc][sampei]
impgc2["ctrl"]<-0
names(impgc2)<-sapply(strsplit(names(impgc2),"_"),"[[",1)
#kdsok<-kdsok[,sampei]

ddf<-melt(cbind(ctrlsok,kdsok))
ddf$type<-"RBP"
ddf$type[grep(ddf$Var2,pattern="ctrl")]<-"control"
ddf$type<-factor(ddf$type,levels = c("RBP","control"))
ddf$type2<-sapply(strsplit(as.character(ddf$Var2),"_"),"[[",1)
ddf$GCimp<-impgc2[as.character(ddf$type2)]
ddf$region<-"end"
ddf<-rbind(ddfs,ddf)
ddf$region<-factor(ddf$region,levels=c("start","end"))
ddf$typetx<-paste("degraded",ex_deg,sep="_")
ddfdeg<-ddf


prof_start<-do.call(lapply(resprofs,function(x){x$profiles[[ex_unc]][500:750,]}),what = cbind)
prof_end<-do.call(lapply(resprofs,function(x){x$profilesend[[ex_unc]][700:950,]}),what = cbind)


ctrls<-prof_start[,colnames(prof_start)=="ctrl"]
kds<-prof_start[,colnames(prof_start)=="kd"]
dups<-duplicated(t(ctrls))
ctrls<-ctrls[,!dups]


colnames(ctrls)<-paste(colnames(ctrls),1:length(colnames(ctrls)),sep = "_")
colnames(kds)<-rbps

impgc<-setNames(df$importance,df$study)[sapply(strsplit(names(resprofs),"_"),function(x){paste(x[-4],collapse="_")})]
ordgc<-order(impgc,decreasing = F)

kdsok<-kds[,ordgc]
set.seed(665)
ctrlsok<-ctrls
#ctrlsok<-ctrls[,sample(1:length(ctrls[1,]),size = 10,replace = F)]

set.seed(666)
impgc2<-impgc
#sampei<-sort(sample(1:length(kdsok[1,]),size = 40,replace = F))
#impgc2<-impgc[ordgc][sampei]
impgc2["ctrl"]<-0
names(impgc2)<-sapply(strsplit(names(impgc2),"_"),"[[",1)
#kdsok<-kdsok[,sampei]

ddf<-melt(cbind(ctrlsok,kdsok))
ddf$type<-"RBP"
ddf$type[grep(ddf$Var2,pattern="ctrl")]<-"control"
ddf$type<-factor(ddf$type,levels = c("RBP","control"))
ddf$type2<-sapply(strsplit(as.character(ddf$Var2),"_"),"[[",1)
ddf$GCimp<-impgc2[as.character(ddf$type2)]
ddf$region<-"start"
ddfs<-ddf

ctrls<-prof_end[,colnames(prof_end)=="ctrl"]
kds<-prof_end[,colnames(prof_end)=="kd"]
dups<-duplicated(t(ctrls))
ctrls<-ctrls[,!dups]


colnames(ctrls)<-paste(colnames(ctrls),1:length(colnames(ctrls)),sep = "_")
colnames(kds)<-rbps

impgc<-setNames(df$importance,df$study)[sapply(strsplit(names(resprofs),"_"),function(x){paste(x[-4],collapse="_")})]
ordgc<-order(impgc,decreasing = F)

kdsok<-kds[,ordgc]
set.seed(665)
ctrlsok<-ctrls
#ctrlsok<-ctrls[,sample(1:length(ctrls[1,]),size = 10,replace = F)]

set.seed(666)
impgc2<-impgc
#sampei<-sort(sample(1:length(kdsok[1,]),size = 40,replace = F))
#impgc2<-impgc[ordgc][sampei]
impgc2["ctrl"]<-0
names(impgc2)<-sapply(strsplit(names(impgc2),"_"),"[[",1)
#kdsok<-kdsok[,sampei]

ddf<-melt(cbind(ctrlsok,kdsok))
ddf$type<-"RBP"
ddf$type[grep(ddf$Var2,pattern="ctrl")]<-"control"
ddf$type<-factor(ddf$type,levels = c("RBP","control"))
ddf$type2<-sapply(strsplit(as.character(ddf$Var2),"_"),"[[",1)
ddf$GCimp<-impgc2[as.character(ddf$type2)]
ddf$region<-"end"
ddf<-rbind(ddfs,ddf)
ddf$region<-factor(ddf$region,levels=c("start","end"))
ddf$typetx<-paste("unchanging",ex_unc,sep="_")
ddfunc<-ddf
ddf<-rbind(ddfst,ddfdeg,ddfunc)

ddf<-ddf[grep(ddf$Var2,pattern = "ctrl",invert = T),]
#collsi<-alpha(c(rep("grey",10),colorRampPalette(c("dark red","grey","dark blue"))(40)),.9)
ddf$importance_GC<-setNames(pvsok$importance_GCpct_cds,pvsok$study)[ddf$Var2]
#aggregate when you got the time!
ddf$typetx<-factor(ddf$typetx,levels=rev(c(paste("degraded",ex_deg,sep="_"),paste("unchanging",ex_unc,sep="_"),paste("stabilized",ex_stab,sep="_"))))

impgc2<-impgc
names(impgc2)<-sapply(strsplit(names(impgc),"_"),"[[",1)
ddf$importance_GC<-impgc2[as.character(ddf$Var2)]
ddf$importance_GCok<-cut(ddf$importance_GC,breaks = quantile(ddf$importance_GC,prob=c(seq(0,1,length.out = 5)),na.rm=T),include.lowest = T,right = T)
ddf$importance_GCok<-gsub(gsub(gsub(gsub(as.character(ddf$importance_GCok),pattern = "[(]",replacement = ""),pattern = "\\[",replacement = ""),pattern = "\\]",replacement = ""),pattern = ",",replacement = "-")


collsi<-alpha(colorRampPalette(c("dark grey","dark blue"))(length(table(as.character(ddf$Var2)))),.8)
collsi2<-alpha(colorRampPalette(c("dark grey","dark blue"))(length(table(ddf$importance_GCok))),.8)


a<-ggplot(ddf,aes(x=Var1,y=value,group=Var2,color=Var2)) 
a<-a + theme_classic() + geom_line() +
  #scale_size_continuous(breaks = c(30,150,300),labels = c(30,150,300),limits = c(20,1500))+
  ylab("coverage") +
  xlab("position") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
  scale_color_manual(values = collsi,"RBP_type") + theme(legend.position = "none")


impgc2<-impgc
names(impgc2)<-sapply(strsplit(names(impgc),"_"),"[[",1)

a<-ggplot(ddf,aes(x=Var1,y=value,group=importance_GCok,color=importance_GCok)) 
a<-a + theme_classic() + stat_summary(geom="line",size=2) +
  #scale_size_continuous(breaks = c(30,150,300),labels = c(30,150,300),limits = c(20,1500))+
  ylab("coverage") +
  xlab("position") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
  scale_color_manual(values = collsi2,"importance_GC") 

#pdf(file = "Figures/Fig5_5cov_examplesENCall_agg.pdf",width = 15,height = 10)
#a + facet_wrap(typetx~region,nrow = 3,scales="free") 
#dev.off()

ddf2<-ddf[ddf$typetx!=paste("unchanging",ex_unc,sep="_"),]


a<-ggplot(ddf2,aes(x=Var1,y=value,group=Var2,color=Var2)) 
a<-a + theme_classic() + geom_line() +
  #scale_size_continuous(breaks = c(30,150,300),labels = c(30,150,300),limits = c(20,1500))+
  ylab("coverage") +
  xlab("position") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
  scale_color_manual(values = collsi,"RBP_type") + theme(legend.position = "none")

#pdf(file = "Figures/Fig5_5cov_examplesENCnounch.pdf",width = 15,height = 8)
#a + facet_wrap(typetx~region,nrow = 3,scales="free") 
#dev.off()

a<-ggplot(ddf2,aes(x=Var1,y=value,group=importance_GCok,color=importance_GCok)) 
a<-a + theme_classic() + stat_summary(geom="line",size=2) +
  #scale_size_continuous(breaks = c(30,150,300),labels = c(30,150,300),limits = c(20,1500))+
  ylab("coverage") +
  xlab("position") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
  scale_color_manual(values = collsi2,"importance_GC") 

#pdf(file = "Figures/Fig5_5cov_examplesENCnounch_agg.pdf",width = 15,height = 8)
#a + facet_wrap(typetx~region,nrow = 3,scales="free") 
#dev.off()

ddf3<-ddf[ddf$typetx!=paste("unchanging",ex_unc,sep="_") & ddf$region=="start",]

collsi<-alpha(colorRampPalette(c("dark grey","dark blue"))(length(unique(ddf3$Var2))),.6)
ddf3$Var2<-factor(ddf3$Var2,levels=unique(ddf3$Var2[order(ddf3$GCimp,decreasing = F)]))
a<-ggplot(ddf3,aes(x=Var1,y=value,group=GCimp,color=GCimp)) 
a<-a + theme_classic() + geom_line() +
  #scale_size_continuous(breaks = c(30,150,300),labels = c(30,150,300),limits = c(20,1500))+
  ylab("coverage") +
  xlab("position") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))+ 
  scale_color_gradient(low = "dark grey",high = "dark blue","Importance\nGCcds")

#pdf(file = "Figures/Fig5_5cov_examplesENCnounch_stonly.pdf",width = 10,height = 8)
#a + facet_wrap(typetx~region,nrow = 3,scales="free") 
#dev.off()


ddf4<-ddf3

ddf4$importance_GCpct_cds<-cut(ddf4$importance_GC,breaks = quantile(pvsok$importance_GCcds,prob=c(seq(0,1,length.out = 5))),include.lowest = T,right = T)

a<-ggplot(ddf4,aes(x=Var1,y=value,group=importance_GCpct_cds,color=importance_GCpct_cds)) 
a<-a + theme_classic() + stat_summary(geom="line",size=2) +
    #scale_size_continuous(breaks = c(30,150,300),labels = c(30,150,300),limits = c(20,1500))+
    ylab("coverage") +
    xlab("position") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
    scale_color_manual(values = collsi2,"importance_GCcds") 

pdf(file = "Figures/Fig6D.pdf",width = 10,height = 8)
a + facet_wrap(typetx~.,nrow = 3,scales="free") 
dev.off()





rf_tab<-read.table("data/RF_table_48h.tsv",header = T,sep = "\t",stringsAsFactors = F)
corse<-c()
for(i in colnames(rf_tab)[3:length(colnames(rf_tab))]){
    corse<-c(corse,cor.test(rf_tab[,2],rf_tab[,i])$estimate)
}
names(corse)<-colnames(rf_tab)[3:length(colnames(rf_tab))]
rf_tab2<-rf_tab
rf_tab2$gene_id<-NULL

limo<-lm(delta_TE~.,data = rf_tab2)

#rf_tab3<-rf_tab2[,colnames(rf_tab2)%in%names(corse[abs(corse)>.15])]
pdf(file = "Figures/Supplementary_Figure_3.pdf",width = 10,height = 10)
corrplot(cor(rf_tab2),method = "color",diag = FALSE,order="hclust",type="upper"
         ,col=colorRampPalette(c("blue","cornflowerblue","white","orange","red"))(10),tl.col = "black",tl.cex = .5)
dev.off()

# make QC plots for Suppl 1
qc1<-get(load("data/QCmetrics.RData"))
qc2<-as.data.frame(get(load("data/QCmetrics2.RData")))
qc2$dataset<-rownames(qc2)
qc2<-qc2[,c(5,1,2,3,4)]
write.table(qc2,file = "Supplementary_Table_1.tsv",sep="\t",quote = F,row.names = F)

qc2<-get(load("data/QCmetrics2.RData"))
qc2<-qc2[,c(4,3,2,1)]
qc2[,4]<-qc2[,4]-qc2[,3]
qc2[,3]<-qc2[,3]-qc2[,2]
qc2[,2]<-qc2[,2]-qc2[,1]

qc2<-melt(qc2)
qc2$seqtype<-"RNA-seq"
qc2$seqtype[grep(qc2$Var1,pattern = "Ribo")]<-"Ribo-seq"
colsi<-rev(colorRampPalette(c("light grey","black"))(4))
g2 <- ggplot(data = qc2, aes(x=Var1, y=value, fill=Var2)) +theme_classic() +
    geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
    labs(x="", y="read count\n")  + theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    scale_fill_manual(values = colsi,"")
g2<-g2 + facet_wrap(~seqtype,nrow = 1,scales = "free")


qc2reg<-do.call(lapply(qc1,function(x){x[["reads_region"]]}),what = rbind)
qc2reg<-qc2reg[c(11,12,1:10),]
qc2reg<-t(apply(qc2reg,1,function(x){x/sum(x)}))*100
qc2reg<-melt(qc2reg)
#colsi<-rev(colorRampPalette(c("light grey","black"))(4))
g2reg <- ggplot(data = qc2reg, aes(x=Var2, y=value, fill=Var2)) +theme_classic() +
    stat_summary(geom = "bar") +
    labs(x="", y="read count %")  + theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    scale_fill_viridis(discrete=T) + theme(legend.position = "none")


qc2rld<-do.call(lapply(qc1,function(x){
    x<-x[["rld"]]
    xnew<-rep(0,22)
    rlds<-as.numeric(gsub(names(x),pattern = "reads_",replacement = ""))
    x<-x[rlds>17]
    rlds<-as.numeric(gsub(names(x),pattern = "reads_",replacement = ""))
    x<-x[rlds<41]
    rlds<-as.numeric(gsub(names(x),pattern = "reads_",replacement = ""))
    x/sum(x)*100
}),what = rbind)
qc2rld<-qc2rld[c(11,12,1:10),]
#qc2rld<-t(apply(qc2rld,1,function(x){x/sum(x)}))*100
qc2rld<-melt(qc2rld)

qc2rld$Var2<-as.numeric(gsub(as.character(qc2rld$Var2),pattern = "reads_",replacement = ""))
#colsi<-rev(colorRampPalette(c("light grey","black"))(4))
g2rld <- ggplot(data = qc2rld, aes(x=Var2, y=value)) +theme_classic() +
    stat_summary() + stat_summary(geom = "line") +
    labs(x="footprint length", y="read count %")  + theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  
#scale_fill_viridis(discrete=T) + theme(legend.position = "none")

#pdf(file = "Figures/Supp_Fig12.pdf",width = 15,height = 12)

#dev.off()

qc2pro<-do.call(lapply(qc1,function(x){
    x<-x[["profile"]]
    unname(x/sum(x)*100)
}),what = rbind)
qc2pro<-qc2pro[c(11,12,1:10),]
#qc2rld<-t(apply(qc2rld,1,function(x){x/sum(x)}))*100
qc2pro<-melt(qc2pro)
res_labels <- c("TSS", "start\ncodon", "stop\ncodon", "TES")
plc<-names(qc1$Ribo_04h_rep1$profile)
res_breaks<-c(grep(plc,pattern = "_1$"),length(plc))

#qc2rld$Var2<-as.numeric(gsub(as.character(qc2rld$Var2),pattern = "reads_",replacement = ""))
#colsi<-rev(colorRampPalette(c("light grey","black"))(4))
g2pro <- ggplot(data = qc2pro, aes(x=Var2, y=value)) +theme_classic() +
    stat_summary()  + scale_x_continuous(breaks=res_breaks, labels=res_labels) +
    labs(x="", y="read count %")  + theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
#scale_fill_viridis(discrete=T) + theme(legend.position = "none")

#pdf(file = "Figures/Supp_Fig12.pdf",width = 15,height = 12)

#dev.off()


pdf(file = "Figures/Supplementary_Figure_1.pdf",width = 18,height = 5)
plot_grid(g2reg,g2rld,g2pro,align = "h",nrow = 1,rel_widths = c(2,2,3))
dev.off()

#cell cycle Albert


df <- read.csv("data/table_031122_ok.tsv",sep = "\t",header = F)
colnames(df)<-c("name","S","G2","G1")
df<-rbind(df[1:3,],df)
df[1:3,1]<-"0hr_Auxin"
df<-melt(df)
df$Treatment<-"DMSO"
df$Treatment[grep(df$name,pattern = "Auxin")]<-"Auxin"
df$Time<-gsub(sapply(strsplit(df$name,"_"),"[[",1),pattern = "hr",replacement = "")
df<-df[order(df$Time),]
df<-rbind(df[1:6,],df)
df[1:6,"Treatment"]<-"Auxin"
df$Time<-as.numeric(df$Time)
df$phtr<-interaction(df$variable,df$Treatment)
df$Phase<-factor(df$variable,levels=c("G1","S","G2"))
ggcc<-ggplot(df, aes(x=Time, y=value, group=phtr, color=Phase,alpha=Treatment)) + 
    stat_summary(geom="line") +
    stat_summary() +
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    ylab("% of total cells") +
    xlab("Hours") +  scale_color_manual(values = c("cornflowerblue","blue","dark blue"),"Cell cycle phase") +
    scale_alpha_manual(values=c("DMSO"=.3,"Auxin"=1)) + scale_x_continuous(breaks = c(0,4,8,16,24,48),labels = as.character(c(0,4,8,16,24,48)))

pdf(file = "Figures/Fig4F.pdf",width = 9,height = 4)
ggcc 
dev.off()

ddx<-ddx_hum


load("data/rleall_stps.RData")
load("data/rleall_stps1.RData")
load("data/rleall_stps2.RData")

reslm<-lapply(runlsok,function(x){
  x<-data.frame(x$pospointone)
  x$time<-c(0,4,8,16,24,48)
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time,data = x)$coefficients[2])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time,data = x)$coefficients[2])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  rescoef
  
})
reslm1<-lapply(runlsok1,function(x){
  x<-data.frame(x$pospointone)
  x$time<-c(0,4,8,16,24,48)
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time,data = x)$coefficients[2])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time,data = x)$coefficients[2])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  rescoef
  
})
reslm2<-lapply(runlsok2,function(x){
  x<-data.frame(x$pospointone)
  x$time<-c(0,4,8,16,24,48)
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time,data = x)$coefficients[2])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time,data = x)$coefficients[2])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  rescoef
  
})

meancovst<-lapply(runlsok, function(x){
  topad<-300-(dim(x$profilesavg)[1]-(x$startp_rna+500))
  if(topad>0){x$profilesavg<-rbind(x$profilesavg,matrix(0,nrow = topad,ncol=ncol(x$profilesavg)))}
  #x$profilesavg<-apply(x$profilesavg,2,function(x){x/max(x,na.rm = T)})
  y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
  
  y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
  x<-data.frame(y)
  x$time<-c(0,4,8,16,24,48)
  
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(x[1,1]==0){x[1,1]<-x[2,1]}
  if(x[1,2]==0){x[1,2]<-x[2,2]}
  x[,1]<-log2(x[,1]/x[1,1])
  x[,2]<-log2(x[,2]/x[1,2])
  if(length(which(is.infinite(x[,1])))>0){x[which(is.infinite(x[,1])),1]<-0}
  if(length(which(is.infinite(x[,2])))>0){x[which(is.infinite(x[,2])),2]<-0}
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time -1,data = x)$coefficients[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time -1,data = x)$coefficients[1])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  resl<-list(y,rescoef)
  names(resl)<-c("avg_cov","lmres")
  resl
  
})
meancovst1<-lapply(runlsok1, function(x){
  topad<-300-(dim(x$profilesavg)[1]-(x$startp_rna+500))
  if(topad>0){x$profilesavg<-rbind(x$profilesavg,matrix(0,nrow = topad,ncol=ncol(x$profilesavg)))}
  #x$profilesavg<-apply(x$profilesavg,2,function(x){x/max(x,na.rm = T)})
  y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
  
  y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
  x<-data.frame(y)
  x$time<-c(0,4,8,16,24,48)
  
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(x[1,1]==0){x[1,1]<-x[2,1]}
  if(x[1,2]==0){x[1,2]<-x[2,2]}
  x[,1]<-log2(x[,1]/x[1,1])
  x[,2]<-log2(x[,2]/x[1,2])
  if(length(which(is.infinite(x[,1])))>0){x[which(is.infinite(x[,1])),1]<-0}
  if(length(which(is.infinite(x[,2])))>0){x[which(is.infinite(x[,2])),2]<-0}
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time -1,data = x)$coefficients[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time -1,data = x)$coefficients[1])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  resl<-list(y,rescoef)
  names(resl)<-c("avg_cov","lmres")
  resl
  
})
meancovst2<-lapply(runlsok2, function(x){
  topad<-300-(dim(x$profilesavg)[1]-(x$startp_rna+500))
  if(topad>0){x$profilesavg<-rbind(x$profilesavg,matrix(0,nrow = topad,ncol=ncol(x$profilesavg)))}
  #x$profilesavg<-apply(x$profilesavg,2,function(x){x/max(x,na.rm = T)})
  y<-x$profilesavg[(x$startp_rna+500-150):(x$startp_rna+500+149),]
  
  y<-matrix(unname(colMeans(y,na.rm = T)),ncol = 2)
  x<-data.frame(y)
  x$time<-c(0,4,8,16,24,48)
  
  colnames(x)<-c("RNA","Ribo","time")
  rescoef<-c(NA,NA,NA,NA)
  if(x[1,1]==0){x[1,1]<-x[2,1]}
  if(x[1,2]==0){x[1,2]<-x[2,2]}
  x[,1]<-log2(x[,1]/x[1,1])
  x[,2]<-log2(x[,2]/x[1,2])
  if(length(which(is.infinite(x[,1])))>0){x[which(is.infinite(x[,1])),1]<-0}
  if(length(which(is.infinite(x[,2])))>0){x[which(is.infinite(x[,2])),2]<-0}
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[1]<-as.numeric(lm(formula = RNA ~ time -1,data = x)$coefficients[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[2]<-as.numeric(lm(formula = Ribo ~ time -1,data = x)$coefficients[1])
  }
  
  if(sum(!is.na(x[,1]))>1){
    rescoef[3]<-as.numeric(anova(lm(formula = RNA ~ time,data = x))$`Pr(>F)`[1])
  }
  if(sum(!is.na(x[,2]))>1){
    rescoef[4]<-as.numeric(anova(lm(formula = Ribo ~ time,data = x))$`Pr(>F)`[1])
  }
  
  resl<-list(y,rescoef)
  names(resl)<-c("avg_cov","lmres")
  resl
  
})


coefsRNA<-sapply(reslm,"[[",1)
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("degraded", "unchanging") )
colsi2<-c("dark grey","blue","red")
pl<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(coverage start position)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = -1.5+1:3,method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


#pdf(file = "Figures/Fig5_lmslope.pdf",width = 10,height = 5)
pl<- pl + coord_cartesian(ylim = c(-6,4))
#dev.off()
coefsRNA<-sapply(meancovst,function(x){x$lmres[1]})
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
#my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("unchanging", "degraded") )
colsi2<-c("dark grey","blue","red")
ppl<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(log2FC coverage)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 0.02+c(.02,.03,.04),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))



coefsRNA<-sapply(reslm1,"[[",1)
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("degraded", "unchanging") )
colsi2<-c("dark grey","blue","red")
pl1<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(coverage start position)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = -1.5+1:3,method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


#pdf(file = "Figures/Fig5_lmslope.pdf",width = 10,height = 5)
pl1<- pl1 + coord_cartesian(ylim = c(-6,4))
#dev.off()

coefsRNA<-sapply(meancovst1,function(x){x$lmres[1]})
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
#my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("unchanging", "degraded") )
colsi2<-c("dark grey","blue","red")
ppl1<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(log2FC coverage)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 0.02+c(.02,.03,.04),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))



coefsRNA<-sapply(reslm2,"[[",1)
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("degraded", "unchanging") )
colsi2<-c("dark grey","blue","red")
pl2<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(coverage start position)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = -1.5+1:3,method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


#pdf(file = "Figures/Fig5_lmslope.pdf",width = 10,height = 5)
pl2<- pl2 + coord_cartesian(ylim = c(-6,4))
#dev.off()

coefsRNA<-sapply(meancovst2,function(x){x$lmres[1]})
ginsclass<-ddx$type_toRle[match(names(coefsRNA),ddx$gene_name)]
dfok<-cbind.data.frame(coefsRNA,ginsclass)
colnames(dfok)<-c("lm_slope","tx_class")
dfok<-dfok[!is.na(dfok$tx_class),]
dfok$tx_class<-factor(dfok$tx_class,levels = c("unchanging","stabilized","degraded"))
#my_comparisons <- list( c("unchanging", "stabilized"), c("degraded", "stabilized"),c("unchanging", "degraded") )
colsi2<-c("dark grey","blue","red")
ppl2<-ggplot(dfok,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(log2FC coverage)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = 0.02+c(.02,.03,.04),method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))


pl1$data$dataset=".10"
pl$data$dataset=".15"
pl2$data$dataset=".20"
plall<-rbind(pl1$data,pl$data,pl2$data)

plall<-ggplot(plall,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(coverage start position)") +
  scale_fill_manual(values = colsi2) +
  geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,label.y = -1.5+1:3,method = "wilcox.test",na.rm = T,method.args = list(alternative="g"))


#pdf(file = "Fig5_lmslope.pdf",width = 10,height = 5)
plall<-plall + coord_cartesian(ylim = c(-6,4)) + facet_wrap(~dataset,nrow = 1) 
#dev.off()

ppl1$data$dataset=".10"
ppl$data$dataset=".15"
ppl2$data$dataset=".20"
pplall<-rbind(ppl1$data,ppl$data,ppl2$data)

pplall<-ggplot(pplall,aes(tx_class,lm_slope,fill=tx_class)) + geom_violin(scale = "area",draw_quantiles = .5)+
  xlab("") +
  ylab("linear model slope\n(log2FC coverage)") +
  scale_fill_manual(values = colsi2) +
  #geom_hline(yintercept = 0,linetype=2) +
  theme_classic() +
  geom_hline(yintercept = 0,linetype=2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
  theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
  theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",na.rm = T,method.args = list(alternative="l"))


#pdf(file = "Fig5_lmslope.pdf",width = 10,height = 5)
pplall<-pplall + facet_wrap(~dataset,nrow = 1)
#dev.off()


#pdf(file = "Figures/Fig5_startcov.pdf",width = 10,height = 5)
#dev.off()

pdf(file = "Figures/Supplementary_Figure_9.pdf",width = 10,height = 10)

plot_grid(plall+ theme(legend.position = "none"),pplall+ theme(legend.position = "none"),nrow = 2,align = "hv")

dev.off()




#RF tables

annotation_file<-"data/gencode.v32.annotation.gtf.gz_Rannot"

all_vars=T

load_annotation(annotation_file)

inspres<-readRDS("data/INSPEcT_res")
gns<-GTF_annotation$trann$gene_name[match(rownames(inspres@ratePvals),GTF_annotation$trann$gene_id)]
#plotGene(inspres, ix=featureNames(inspres)[gns=="TLE4"],priors = F,constantModel = T)
summary_multiDE<-get(load("data/degron_multiDE_results_summary.RData"))

df<-summary_multiDE$uniq$`48h`

red_ex <- unlist(GTF_annotation$exons_txs[names(GTF_annotation$cds_txs)])
red_ex$gene_id<-GTF_annotation$trann$gene_id[match(names(red_ex),GTF_annotation$trann$transcript_id)]
red_ex<-reduce(split(red_ex,red_ex$gene_id))

gcexcds<-letterFrequency(extractTranscriptSeqs(x=genome_seq,red_ex),letters = "GC",as.prob = T)
gcexcds<-setNames(gcexcds[,1],names(red_ex))
df$GCpct_ex<-gcexcds[df$gene_id]
#newres_DE_DEX$annotation_files<-list_annot
codop<-read.table("data/13059_2020_2251_MOESM6_ESM.csv",sep = ",",stringsAsFactors = F,header = T,skip=1)
codop<-codop[codop$specie=="human",]

gnsi<-sapply(strsplit(df$gene_id,"[.]"),"[[",1)
mtt<-match(gnsi,codop$gene_id)
ok<-!is.na(mtt)
df$codonopt1<-NA
df$codonopt2<-NA
df$codonopt1[ok]<-codop$PLS1[mtt[ok]]
df$codonopt2[ok]<-codop$PLS2[mtt[ok]]

df$TPM_RNA<-as.numeric(df$TPM_Intron)
df$TPM_Ribo<-as.numeric(df$TPM_Intron)
df$TPM_Intron<-as.numeric(df$TPM_Intron)

df[,colnames(df)[grep(colnames(df),pattern = "len")]]<-log2(df[,colnames(df)[grep(colnames(df),pattern = "len")]]+1)
df[,colnames(df)[grep(colnames(df),pattern = "TPM")]]<-log2(df[,colnames(df)[grep(colnames(df),pattern = "TPM")]])

#df$scores_top<-as.numeric(toupper(df$gene_name)%in%topgns[,1])

dfok<-df[,colnames(df)[grep(colnames(df),pattern = "TPM_RNA|len|pct|^RNA_log2FC|base_|pos|dens|score|delta|n_exons|codonfr|codonopt")]]
dfok<-dfok[,grep(colnames(dfok),pattern = "posdex",invert = T)]
dfok<-dfok[,grep(colnames(dfok),pattern = "binsdex",invert = T)]

rownames(dfok)<-df$gene_id
dfok[which(dfok==Inf,arr.ind = T)]<-NA
dfok[which(dfok==-Inf,arr.ind = T)]<-NA
yeah<-which(dfok$TPM_RNA>3)
if(!length(yeah)>1000){
  dfok$TPM_RNA<-df$TPM_Ribo
  yeah<-which(dfok$TPM_RNA>3)
}

dfok<-dfok[yeah,]

#dfok<-dfok[,which(!(grepl(colnames(dfok),pattern = "pos") & grepl(colnames(dfok),pattern = "FC")))]
#put those NAs to zero and check base during training

correctNA<-dfok[,grep(colnames(dfok),pattern = "base|delta")]
correctNA[is.na(correctNA)]<-0
dfok[,grep(colnames(dfok),pattern = "base|delta")]<-correctNA


totest<-colnames(dfok)[grep(colnames(dfok),pattern = "log2|delta")]
totestall<-totest
if(!all_vars){totest<-c("RNA_log2FC","delta_TE","delta_intrex")}
checkok<-rep(T,length(totest))
for(toti in 1:length(totest)){
  checkok[toti]<-length(unique(dfok[,totest[toti]]))>10
}


xyribo<-setNames(paste("RiboRNA",df$xy_RiboRNA,sep = "_"),df$gene_id)[rownames(dfok)]
xyintron<-setNames(paste("IntronExon",df$xy_IntronExon,sep = "_"),df$gene_id)[rownames(dfok)]

distxyribo<-setNames(df$mindist_RiboRNA,df$gene_id)[rownames(dfok)]
min_xy<-sort(distxyribo,decreasing = T)[as.integer(dim(dfok)[1]/5)]
#Think about ways to improve the xy thing, but probably better not to mess with it,
#needs more brain!
#xyribo[abs(distxyribo)>min_xy]<-"other"

xyribo<-c(split(dfok,xyribo),split(dfok,xyintron))
xyribo[["all"]]<-dfok
xyribo<-xyribo[grep(names(xyribo),pattern = "_NA",invert = T)]
set.seed(666)
sizz<-sample(1:1000,replace = F,size = length(totest))
list_resrff_all<-list()
nam="all"
list_resrff<-list()
dfok<-xyribo[[nam]]
i="delta_TE"
dfok2<-dfok[complete.cases(dfok[,i]),]
if(grepl(i,pattern = "delta")){
  bss<-gsub(i,pattern = "delta",replacement = "base")
  dfok2<-dfok2[dfok2[,bss]!=0,]
}
if(dim(dfok2)[1]<100){next}
dfok2$to_test<-dfok2[,i]

dfok2<-dfok2[,colnames(dfok2)!=i]
onedelta=T
if(as.logical(as.character(onedelta))){
  dfok2<-dfok2[,!colnames(dfok2)%in%totestall]
}


dfok2<-na.roughfix(dfok2)
dfok2$delta_TE<-dfok2$to_test
dfok2$to_test<-NULL
dfok2$gene_id<-rownames(dfok2)
colni<-sort(colnames(dfok2))
colni<-colni[-which(colni%in%c("gene_id","delta_TE"))]
dfok2<-dfok2[,c("gene_id","delta_TE",colni)]
write.table(dfok2,file = "Tables/Supplementary_Table_2.tsv",row.names = F,quote = F,sep="\t")


#mouse table

annotation_file<-"data/gencode.vM22.annotation.gtf_Rannot"

all_vars=T

load_annotation(annotation_file)


summary_multiDE<-get(load("data/mouseAllnoHet_multiDE_results_summary.RData"))

df<-summary_multiDE$uniq$cKO

red_ex <- unlist(GTF_annotation$exons_txs[names(GTF_annotation$cds_txs)])
red_ex$gene_id<-GTF_annotation$trann$gene_id[match(names(red_ex),GTF_annotation$trann$transcript_id)]
red_ex<-reduce(split(red_ex,red_ex$gene_id))

gcexcds<-letterFrequency(extractTranscriptSeqs(x=genome_seq,red_ex),letters = "GC",as.prob = T)
gcexcds<-setNames(gcexcds[,1],names(red_ex))
df$GCpct_ex<-gcexcds[df$gene_id]
#newres_DE_DEX$annotation_files<-list_annot
codop<-read.table("data/13059_2020_2251_MOESM6_ESM.csv",sep = ",",stringsAsFactors = F,header = T,skip=1)
codop<-codop[codop$specie=="mouse",]

gnsi<-sapply(strsplit(df$gene_id,"[.]"),"[[",1)
mtt<-match(gnsi,codop$gene_id)
ok<-!is.na(mtt)
df$codonopt1<-NA
df$codonopt2<-NA
df$codonopt1[ok]<-codop$PLS1[mtt[ok]]
df$codonopt2[ok]<-codop$PLS2[mtt[ok]]

df$TPM_RNA<-as.numeric(df$TPM_Intron)
df$TPM_Ribo<-as.numeric(df$TPM_Intron)
df$TPM_Intron<-as.numeric(df$TPM_Intron)

df[,colnames(df)[grep(colnames(df),pattern = "len")]]<-log2(df[,colnames(df)[grep(colnames(df),pattern = "len")]]+1)
df[,colnames(df)[grep(colnames(df),pattern = "TPM")]]<-log2(df[,colnames(df)[grep(colnames(df),pattern = "TPM")]])

#df$scores_top<-as.numeric(toupper(df$gene_name)%in%topgns[,1])

dfok<-df[,colnames(df)[grep(colnames(df),pattern = "TPM_RNA|len|pct|^RNA_log2FC|base_|pos|dens|score|delta|n_exons|codonfr|codonopt")]]
dfok<-dfok[,grep(colnames(dfok),pattern = "posdex",invert = T)]
dfok<-dfok[,grep(colnames(dfok),pattern = "binsdex",invert = T)]

rownames(dfok)<-df$gene_id
dfok[which(dfok==Inf,arr.ind = T)]<-NA
dfok[which(dfok==-Inf,arr.ind = T)]<-NA
yeah<-which(dfok$TPM_RNA>3)
if(!length(yeah)>1000){
  dfok$TPM_RNA<-df$TPM_Ribo
  yeah<-which(dfok$TPM_RNA>3)
}

dfok<-dfok[yeah,]

#dfok<-dfok[,which(!(grepl(colnames(dfok),pattern = "pos") & grepl(colnames(dfok),pattern = "FC")))]
#put those NAs to zero and check base during training

correctNA<-dfok[,grep(colnames(dfok),pattern = "base|delta")]
correctNA[is.na(correctNA)]<-0
dfok[,grep(colnames(dfok),pattern = "base|delta")]<-correctNA


totest<-colnames(dfok)[grep(colnames(dfok),pattern = "log2|delta")]
totestall<-totest
if(!all_vars){totest<-c("RNA_log2FC","delta_TE","delta_intrex")}
checkok<-rep(T,length(totest))
for(toti in 1:length(totest)){
  checkok[toti]<-length(unique(dfok[,totest[toti]]))>10
}


xyribo<-setNames(paste("RiboRNA",df$xy_RiboRNA,sep = "_"),df$gene_id)[rownames(dfok)]
xyintron<-setNames(paste("IntronExon",df$xy_IntronExon,sep = "_"),df$gene_id)[rownames(dfok)]

distxyribo<-setNames(df$mindist_RiboRNA,df$gene_id)[rownames(dfok)]
min_xy<-sort(distxyribo,decreasing = T)[as.integer(dim(dfok)[1]/5)]
#Think about ways to improve the xy thing, but probably better not to mess with it,
#needs more brain!
#xyribo[abs(distxyribo)>min_xy]<-"other"

xyribo<-c(split(dfok,xyribo),split(dfok,xyintron))
xyribo[["all"]]<-dfok
xyribo<-xyribo[grep(names(xyribo),pattern = "_NA",invert = T)]
set.seed(666)
sizz<-sample(1:1000,replace = F,size = length(totest))
list_resrff_all<-list()
nam="all"
list_resrff<-list()
dfok<-xyribo[[nam]]
i="delta_TE"
dfok2<-dfok[complete.cases(dfok[,i]),]
if(grepl(i,pattern = "delta")){
  bss<-gsub(i,pattern = "delta",replacement = "base")
  dfok2<-dfok2[dfok2[,bss]!=0,]
}
if(dim(dfok2)[1]<100){next}
dfok2$to_test<-dfok2[,i]

dfok2<-dfok2[,colnames(dfok2)!=i]
onedelta=T
if(as.logical(as.character(onedelta))){
  dfok2<-dfok2[,!colnames(dfok2)%in%totestall]
}


dfok2<-na.roughfix(dfok2)
dfok2$delta_TE<-dfok2$to_test
dfok2$to_test<-NULL
dfok2$gene_id<-rownames(dfok2)
colni<-sort(colnames(dfok2))
colni<-colni[-which(colni%in%c("gene_id","delta_TE"))]
dfok2<-dfok2[,c("gene_id","delta_TE",colni)]
write.table(dfok2,file = "Tables/Supplementary_Table_4.tsv",row.names = F,quote = F,sep="\t")

#last part, do suppl3

tibook<-read.table("data/all_encode_files.txt",header = T,sep = "\t",stringsAsFactors = F)

write.table(tibook,file = "Tables/Supplementary_Table_3.tsv",row.names = F,quote = F,sep="\t")


