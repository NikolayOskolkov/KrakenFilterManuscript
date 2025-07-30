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


#sessionInfo()
#R version 4.2.3 (2023-03-15)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.6 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=sv_SE.UTF-8        LC_COLLATE=en_US.UTF-8    
#[5] LC_MONETARY=sv_SE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=sv_SE.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=sv_SE.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggplot2_3.4.3      viridis_0.6.3      viridisLite_0.4.2  matrixStats_0.63.0 pheatmap_1.0.12   

#loaded via a namespace (and not attached):
#  [1] rstudioapi_0.14    magrittr_2.0.3     tidyselect_1.2.0   munsell_0.5.0      colorspace_2.1-0   R6_2.5.1           rlang_1.1.6       
#[8] fansi_1.0.5        dplyr_1.1.2        tools_4.2.3        grid_4.2.3         gtable_0.3.4       utf8_1.2.4         cli_3.6.5         
#[15] withr_2.5.0        tibble_3.2.1       lifecycle_1.0.4    gridExtra_2.3      farver_2.1.1       RColorBrewer_1.1-3 vctrs_0.6.5       
#[22] glue_1.6.2         labeling_0.4.2     compiler_4.2.3     pillar_1.9.0       generics_0.1.3     scales_1.3.0       pkgconfig_2.0.3   


