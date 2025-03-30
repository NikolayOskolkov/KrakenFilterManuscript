################################################ GROUND TRUTH MAIN AMETA ANALYSIS ###################################################
library("pheatmap")
gt_path<-"/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/data_and_lists_for_16_species_simulation/"
ground_truth_matrix<-read.delim(paste0(gt_path,"ground_truth_number_of_reads.txt"),row.names=1,header=TRUE,sep="\t")
#rownames(ground_truth_matrix)<-gsub("\\.fa","",rownames(ground_truth_matrix))
ground_truth_matrix<-ground_truth_matrix[order(rownames(ground_truth_matrix)),]
#colnames(ground_truth_matrix)<-paste0("Sample",seq(1,10,1))
pheatmap(ground_truth_matrix,display_numbers=TRUE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Ground truth: sedaDNA samples",number_format="%i")
ground_truth_matrix[1:5,]
colSums(ground_truth_matrix)

ground_truth_matrix_binary<-ground_truth_matrix
ground_truth_matrix_binary[ground_truth_matrix_binary>0]<-1
ground_truth_matrix_binary[ground_truth_matrix_binary<=0]<-0
pheatmap(ground_truth_matrix_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,legend_breaks = c(0,1),
         main="Binary ground truth: sedaDNA samples")
ground_truth_matrix_binary
colSums(ground_truth_matrix_binary)


################################################# KRAKENUNIQ MAIN AMETA ANALYSIS ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")
max_kmers<-10000; step_kmers<-50
max_taxReads<-300; step_taxReads<-10
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
taxReads_vector<-seq(from=0,to=max_taxReads,by=step_taxReads)
sample_vector<-seq(from=1,to=2,by=1)
F1_array<-array(rep(NA, length(sample_vector)*length(taxReads_vector)*length(kmers_vector)),
                c(length(sample_vector), length(taxReads_vector), length(kmers_vector)))
for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  F1_matrix<-matrix(NA,nrow=length(kmers_vector),ncol=length(taxReads_vector))
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  for(j in 1:length(taxReads_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>taxReads_vector[j],]
    #krakenuniq<-krakenuniq[krakenuniq$reads>taxReads_vector[j],]
    for(i in 1:length(kmers_vector))
    {
      krakenuniq<-krakenuniq[krakenuniq$kmers>kmers_vector[i],]
      
      query_list<-krakenuniq$taxName
      true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
      TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
      F1<-(2*TP)/(2*TP+FP+FN)
      
      F1_matrix[i,j]<-F1
      F1_array[s,j,i]<-F1
    }
  }
  colnames(F1_matrix)<-taxReads_vector
  rownames(F1_matrix)<-kmers_vector
  print(head(F1_matrix))
  
  filled.contour(F1_matrix,plot.axes= {
    axis(2,at=as.numeric(colnames(F1_matrix))/max_taxReads,labels=as.numeric(colnames(F1_matrix)))
    axis(1,at=as.numeric(rownames(F1_matrix))/max_kmers,las=2,labels=as.numeric(rownames(F1_matrix)))},
    nlevels=20,color.palette=terrain.colors,main=paste0("Sample",sample_vector[s]))
  
  print(which(F1_matrix==max(F1_matrix), arr.ind=T))
  mtext(paste0("Optimal kmers = ",rownames(which(F1_matrix==max(F1_matrix), arr.ind=T))[1],", F1_max = ",max(F1_matrix)))
}
mean_across_samples<-apply(F1_array,c(2,3),mean)
sd_across_samples<-apply(F1_array,c(2,3),sd)
rownames(mean_across_samples)<-taxReads_vector
colnames(mean_across_samples)<-kmers_vector
rownames(sd_across_samples)<-taxReads_vector
colnames(sd_across_samples)<-kmers_vector

mean_across_samples
max(mean_across_samples)
which(mean_across_samples==max(mean_across_samples), arr.ind=T)

sd_across_samples

mean_across_samples[1,29]
sd_across_samples[1,29]

filled.contour(mean_across_samples,plot.axes= {
  axis(2,at=as.numeric(colnames(mean_across_samples))/max_kmers,labels=as.numeric(colnames(mean_across_samples)))
  axis(1,at=as.numeric(rownames(mean_across_samples))/max_taxReads,las=2,labels=as.numeric(rownames(mean_across_samples)))},
  nlevels=30,color.palette=terrain.colors,main="F1-score average across samples: sedaDNA dataset",xlab="Number of taxReads",
  ylab="Number of unique k-mers")

mtext(paste0("kmers_max = ",colnames(mean_across_samples)[which(mean_across_samples==max(mean_across_samples), arr.ind=T)[,"col"]][1],
             ", taxReads_max = ",rownames(mean_across_samples)[which(mean_across_samples==max(mean_across_samples), arr.ind=T)[,"row"]][1],
             ", F1_max = ",max(mean_across_samples)))




################################################# KRAKENUNIQ FILTER K ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")
par(mfrow=c(2,3))

max_kmers<-10000; step_kmers<-100
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
sample_vector<-seq(from=1,to=2,by=1)

F1_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(kmers_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[krakenuniq$kmers>kmers_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix[i,s]<-F1
  }
  colnames(F1_matrix)<-sample_vector
  rownames(F1_matrix)<-kmers_vector
  print(head(F1_matrix))
}
plot(rowMeans(F1_matrix)~kmers_vector,type="o",xlab="K",
     ylab="F1 of ground truth reconstruction",col="darkblue",ylim=c(0,1),pch=19,main="KrakenUniq K filter")
arrows(x0=kmers_vector, y0=rowMeans(F1_matrix)-rowSds(F1_matrix), x1=kmers_vector, y1=rowMeans(F1_matrix)+rowSds(F1_matrix), 
       code=3, angle=90, length=0.03, col="darkblue")

mtext(paste0("kmers_max = ",names(rowMeans(F1_matrix))[rowMeans(F1_matrix)==max(rowMeans(F1_matrix))][1],
             ", F1_max = ",round(max(rowMeans(F1_matrix)),3)))

rowMeans(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]
rowSds(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]

################################################# KRAKENUNIQ FILTER C ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")

kmers_vector<-seq(from=0,to=0.005,by=0.00001)
sample_vector<-seq(from=1,to=2,by=1)

F1_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(kmers_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[krakenuniq$cov>kmers_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix[i,s]<-F1
  }
  colnames(F1_matrix)<-sample_vector
  rownames(F1_matrix)<-kmers_vector
  print(head(F1_matrix))
}
plot(rowMeans(F1_matrix)~kmers_vector,type="o",xlab="C",
     ylab="F1 of ground truth reconstruction",col="magenta",ylim=c(0,1),pch=19,main="KrakenUniq C filter")
arrows(x0=kmers_vector, y0=rowMeans(F1_matrix)-rowSds(F1_matrix), x1=kmers_vector, y1=rowMeans(F1_matrix)+rowSds(F1_matrix), 
       code=3, angle=90, length=0.03, col="magenta")

mtext(paste0("C_max = ",names(rowMeans(F1_matrix))[rowMeans(F1_matrix)==max(rowMeans(F1_matrix))][1],
             ", F1_max = ",round(max(rowMeans(F1_matrix)),3)))

rowMeans(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]
rowSds(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]



################################################# KRAKENUNIQ FILTER R ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")

kmers_vector<-seq(from=0,to=1000,by=10)
sample_vector<-seq(from=1,to=2,by=1)

F1_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(kmers_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[krakenuniq$reads>kmers_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix[i,s]<-F1
  }
  colnames(F1_matrix)<-sample_vector
  rownames(F1_matrix)<-kmers_vector
  print(head(F1_matrix))
}
plot(rowMeans(F1_matrix)~kmers_vector,type="o",xlab="R",
     ylab="F1 of ground truth reconstruction",col="cyan",ylim=c(0,1),pch=19,main="KrakenUniq R filter")
arrows(x0=kmers_vector, y0=rowMeans(F1_matrix)-rowSds(F1_matrix), x1=kmers_vector, y1=rowMeans(F1_matrix)+rowSds(F1_matrix), 
       code=3, angle=90, length=0.03, col="cyan")

mtext(paste0("R_max = ",names(rowMeans(F1_matrix))[rowMeans(F1_matrix)==max(rowMeans(F1_matrix))][1],
             ", F1_max = ",round(max(rowMeans(F1_matrix)),3)))

rowMeans(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]
rowSds(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]


################################################# KRAKENUNIQ FILTER K/R ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")

ratio_vector<-seq(from=0,to=100,by=0.5)
sample_vector<-seq(from=1,to=2,by=1)

F1_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(ratio_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(ratio_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[(krakenuniq$kmers/krakenuniq$taxReads)>ratio_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix[i,s]<-F1
  }
  colnames(F1_matrix)<-sample_vector
  rownames(F1_matrix)<-ratio_vector
  print(head(F1_matrix))
}
plot(rowMeans(F1_matrix)~ratio_vector,type="o",xlab="K / R",
     ylab="F1 of ground truth reconstruction",col="darkred",ylim=c(0,1),pch=19,main="KrakenUniq K / R filter")
arrows(x0=ratio_vector, y0=rowMeans(F1_matrix)-rowSds(F1_matrix), x1=ratio_vector, y1=rowMeans(F1_matrix)+rowSds(F1_matrix), 
       code=3, angle=90, length=0.03, col="darkred")

mtext(paste0("(K / R)_max = ",names(rowMeans(F1_matrix))[rowMeans(F1_matrix)==max(rowMeans(F1_matrix))][1],
             ", F1_max = ",round(max(rowMeans(F1_matrix)),3)))

mean(as.numeric(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),]))
sd(as.numeric(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),]))


################################################# KRAKENUNIQ FILTER (K/R)*C ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")

ratio_vector<-seq(from=0,to=0.1,by=0.0005)
sample_vector<-seq(from=1,to=2,by=1)

F1_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(ratio_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(ratio_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[((krakenuniq$kmers/krakenuniq$taxReads)*krakenuniq$cov)>ratio_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix[i,s]<-F1
  }
  colnames(F1_matrix)<-sample_vector
  rownames(F1_matrix)<-ratio_vector
  print(head(F1_matrix))
}
plot(rowMeans(F1_matrix)~ratio_vector,type="o",xlab="(K / R) * C",
     ylab="F1 of ground truth reconstruction",col="darkgreen",ylim=c(0,1),pch=19,main="KrakenUniq (K / R) * C filter")
arrows(x0=ratio_vector, y0=rowMeans(F1_matrix)-rowSds(F1_matrix), x1=ratio_vector, y1=rowMeans(F1_matrix)+rowSds(F1_matrix), 
       code=3, angle=90, length=0.03, col="darkgreen")

mtext(paste0("((K / R) * C)_max = ",names(rowMeans(F1_matrix))[rowMeans(F1_matrix)==max(rowMeans(F1_matrix))][1],
             ", F1_max = ",round(max(rowMeans(F1_matrix)),3)))

rowMeans(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]
rowSds(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),])[1]




################################################# KRAKENUNIQ FILTER (K/R)*dexp(C) ##################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
library("matrixStats")

ratio_vector<-seq(from=0,to=100,by=1)
sample_vector<-seq(from=1,to=2,by=1)

F1_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(ratio_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(ratio_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[((krakenuniq$kmers/krakenuniq$taxReads)*(1.3^(18^krakenuniq$cov)))>ratio_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix[i,s]<-F1
  }
  colnames(F1_matrix)<-sample_vector
  rownames(F1_matrix)<-ratio_vector
  print(head(F1_matrix))
}
plot(rowMeans(F1_matrix)~ratio_vector,type="o",xlab="(K / R) * dexp(C)",
     ylab="F1 of ground truth reconstruction",col=j,ylim=c(0,1),pch=19,main="KrakenUniq (K / R) * dexp(C) filter")
arrows(x0=ratio_vector, y0=rowMeans(F1_matrix)-rowSds(F1_matrix), x1=ratio_vector, y1=rowMeans(F1_matrix)+rowSds(F1_matrix), 
       code=3, angle=90, length=0.03, col=j)

mtext(paste0("((K / R) * dexp(C))_max = ",names(rowMeans(F1_matrix))[rowMeans(F1_matrix)==max(rowMeans(F1_matrix))][1],
             ", F1_max = ",round(max(rowMeans(F1_matrix)),3)))

mean(as.numeric(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),]))
sd(as.numeric(F1_matrix[which(rowMeans(F1_matrix)==max(rowMeans(F1_matrix))),]))

################################################ CORRELATION BETWEEN QUALITY FILTERS #################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
sample_vector<-seq(from=1,to=2,by=1)
for(s in 1:2)#length(sample_vector)
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_NT_sample",sample_vector[s]),comment.char="#",
                         check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
  
  cor_mat[[s]]<-cor(subset(krakenuniq,select=c("%","reads","taxReads","kmers","cov","dup")),method="spearman")
}
average_cor <- Reduce("+", cor_mat) / length(cor_mat)
pheatmap(average_cor,display_numbers=TRUE)
