# To run the script and reproduce the figures from the manuscript

# Refining filtering criteria for accurate taxonomic classification of ancient metagenomic data
# Nikolay Oskolkov
# bioRxiv 2025.03.31.646431; doi: https://doi.org/10.1101/2025.03.31.646431 

# please clone the repository https://github.com/NikolayOskolkov/KrakenFilterManuscript



#setwd("/home/nikolay/WABI/Frontiers_in_Microbiology/github/KrakenFilterManuscript/krakenuniq_real_datasets/")
setwd("path_to_cloned_github_repo/krakenuniq_real_datasets/")

file_list_krakenuniq<-c("gok2c_krakenuniq.output_Full_NT",
                        "adycha_krakenuniq.output_Full_NT",
                        "krestovka_krakenuniq.output_Full_NT",
                        "stg001_krakenuniq.output_Full_NT",
                        "ans004_krakenuniq.output_Full_NT",
                        "blank_krakenuniq.output_Full_NT")

gt<-read.delim("krakenuniq_real_samples_ground_truth_matrix_binary.txt",header=TRUE,
               check.names=FALSE,row.names=1,sep="\t")
#colnames(gt)<-gsub("_krakenuniq.output_Full_NT","",colnames(gt))
gt

library("pheatmap")
pheatmap(gt,display_numbers=FALSE,fontsize=15,cluster_cols=FALSE,cluster_rows=FALSE,legend_breaks = c(0,1),
         main="Binary ground truth: real samples")

library("matrixStats")

max_kmers<-5000; step_kmers<-100
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
sample_vector<-seq(from=1,to=length(file_list_krakenuniq),by=1)

F1_matrix_ku<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(i in 1:length(file_list_krakenuniq))
{
  krakenuniq<-read.delim(file_list_krakenuniq[i],header=TRUE,comment.char="#",sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  names(krakenuniq)[1]<-"Pers_Reads"
  #krakenuniq<-krakenuniq[krakenuniq$taxReads>50,]
  #krakenuniq<-krakenuniq[krakenuniq$dup<1.1,]
  krakenuniq<-krakenuniq[krakenuniq$cov>0.0015,]
  
  for(j in 1:length(kmers_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$kmers>kmers_vector[j],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(gt)[gt[,file_list_krakenuniq[i]]==1]
    
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix_ku[j,i]<-F1
  }
  colnames(F1_matrix_ku)<-sample_vector
  rownames(F1_matrix_ku)<-kmers_vector
  print(head(F1_matrix_ku))
}
par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(rowMeans(F1_matrix_ku)~kmers_vector,type="o",xlab="Number of unique k-mers",
     ylab="F1 of ground truth reconstruction",col="darkblue",ylim=c(0,0.3),pch=19,
     main="KrakenUniq number of unique k-mers filter: real datasets",
     cex.main=2.2,cex.lab=2.2)
arrows(x0=kmers_vector, y0=rowMeans(F1_matrix_ku)-rowSds(F1_matrix_ku), x1=kmers_vector, y1=rowMeans(F1_matrix_ku)+rowSds(F1_matrix_ku), 
       code=3, angle=90, length=0.03, col="darkblue")





