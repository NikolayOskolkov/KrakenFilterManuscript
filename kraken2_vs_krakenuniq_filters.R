# To run the script and reproduce the figures from the manuscript

# Refining filtering criteria for accurate taxonomic classification of ancient metagenomic data
# Nikolay Oskolkov
# bioRxiv 2025.03.31.646431; doi: https://doi.org/10.1101/2025.03.31.646431 

# please clone the repository https://github.com/NikolayOskolkov/KrakenFilterManuscript


############################################ GROUND TRUTH MAIN AMETA ANALYSIS #############################################

library("pheatmap")

#gt_path<-"/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/data_and_lists_for_main_aMeta_analysis/"
#gt_path<-"/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/data_and_lists_for_pathogen_aMeta_analysis/"
gt_path<-"path_to_cloned_github_repo/krakenuniq_sim_pathogenic_microbial/"

ground_truth_matrix<-read.delim(paste0(gt_path,"ground_truth_number_of_reads.txt"),row.names=1,header=TRUE,sep="\t")
rownames(ground_truth_matrix)<-gsub("\\.fa","",rownames(ground_truth_matrix))
ground_truth_matrix<-ground_truth_matrix[order(rownames(ground_truth_matrix)),]
colnames(ground_truth_matrix)<-paste0("Sample",seq(1,10,1))
pheatmap(ground_truth_matrix,display_numbers=TRUE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Ground truth: microbial species",number_format="%i")
ground_truth_matrix[1:5,]
colSums(ground_truth_matrix)

ground_truth_matrix_binary<-ground_truth_matrix
ground_truth_matrix_binary[ground_truth_matrix_binary>0]<-1
ground_truth_matrix_binary[ground_truth_matrix_binary<=0]<-0
pheatmap(ground_truth_matrix_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         legend_breaks = c(0,1),main="Binary ground truth: microbial species")
ground_truth_matrix_binary[1:5,1:5]
colSums(ground_truth_matrix_binary)


################################################# KRAKENUNIQ FILTER K ####################################################

#ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
#ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_pathogens_unfiltered/"
ku_path<-"path_to_cloned_github_repo/krakenuniq_sim_pathogenic_microbial/"

library("matrixStats")

max_kmers<-20000; step_kmers<-100
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
sample_vector<-seq(from=1,to=10,by=1)

F1_matrix_ku<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",
                         check.names=FALSE,sep="\t")
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
    
    F1_matrix_ku[i,s]<-F1
  }
  colnames(F1_matrix_ku)<-sample_vector
  rownames(F1_matrix_ku)<-kmers_vector
  print(head(F1_matrix_ku))
}
plot(rowMeans(F1_matrix_ku)~kmers_vector,type="o",xlab="K",
     ylab="F1 of ground truth reconstruction",col="darkblue",ylim=c(0,1),pch=19,main="KrakenUniq K filter")
arrows(x0=kmers_vector, y0=rowMeans(F1_matrix_ku)-rowSds(F1_matrix_ku), x1=kmers_vector, y1=rowMeans(F1_matrix_ku)+rowSds(F1_matrix_ku), 
       code=3, angle=90, length=0.03, col="darkblue")

mtext(paste0("kmers_max = ",names(rowMeans(F1_matrix_ku))[rowMeans(F1_matrix_ku)==max(rowMeans(F1_matrix_ku))][1],
             ", F1_max = ",max(rowMeans(F1_matrix_ku))))



################################################# KRAKEN2 FILTER MINIMIZER ##################################################

#k2_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/kraken2_aMeta_simulations_unfiltered/"
#k2_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/kraken2_aMeta_simulations_pathogens_unfiltered/"
k2_path<-"path_to_cloned_github_repo/kraken2_aMeta_simulations_pathogens_unfiltered/"

library("matrixStats")

max_minimizers<-20000; step_minimizers<-100
minimizers_vector<-seq(from=0,to=max_minimizers,by=step_minimizers)
sample_vector<-seq(from=1,to=10,by=1)

F1_matrix_k2<-matrix(NA,ncol=length(sample_vector),nrow=length(minimizers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  kraken2<-read.delim(paste0(k2_path,"sample",sample_vector[s],".trimmed.fastq.gz_kraken2.output"),comment.char="#",
                         check.names=FALSE,header=FALSE,sep="\t")
  colnames(kraken2)<-c("%","reads","taxReads","n_minimizers","n_distinct_minimizers","rank","taxID","taxName")
  kraken2$taxName<-trimws(as.character(kraken2$taxName))
  #kraken2<-kraken2[grepl("S",as.character(kraken2$rank)),]
  kraken2<-kraken2[as.character(kraken2$rank)=="S",]
  
  for(i in 1:length(minimizers_vector))
  {
    kraken2<-kraken2[kraken2$taxReads>0,]
    kraken2<-kraken2[kraken2$n_distinct_minimizers>minimizers_vector[i],]
    #kraken2<-kraken2[kraken2$n_minimizers>minimizers_vector[i],]
    
    query_list<-as.character(kraken2$taxName)
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    
    TP<-sum(query_list%in%true_list); FP<-sum(!query_list%in%true_list); FN<-sum(!true_list%in%query_list)
    F1<-(2*TP)/(2*TP+FP+FN)
    
    F1_matrix_k2[i,s]<-F1
  }
  colnames(F1_matrix_k2)<-sample_vector
  rownames(F1_matrix_k2)<-minimizers_vector
  print(head(F1_matrix_k2))
}
plot(rowMeans(F1_matrix_k2)~minimizers_vector,type="o",xlab="K",
     ylab="F1 of ground truth reconstruction",col="darkblue",ylim=c(0,1),pch=19,main="Kraken2 minimizer filter")
arrows(x0=minimizers_vector, y0=rowMeans(F1_matrix_k2)-rowSds(F1_matrix_k2), x1=minimizers_vector, 
       y1=rowMeans(F1_matrix_k2)+rowSds(F1_matrix_k2), 
       code=3, angle=90, length=0.03, col="darkblue")

mtext(paste0("minimizers_max = ",names(rowMeans(F1_matrix_k2))[rowMeans(F1_matrix_k2)==max(rowMeans(F1_matrix_k2))][1],
             ", F1_max = ",max(rowMeans(F1_matrix_k2))))





#############################################
plot(rowMeans(F1_matrix_ku)~kmers_vector,type="o",xlab="unique kmer / minimizer counts",
     ylab="F1 of ground truth reconstruction",col="darkblue",ylim=c(0,1),pch=19,
     main="KrakenUniq / Kraken2 K filter: pathogenic microbial dataset")
lines(rowMeans(F1_matrix_k2)~minimizers_vector,type="o",col="orange",pch=19)
legend("topright", inset=.02, c("KrakenUniq","Kraken2"), fill=c("darkblue","orange"), cex=1.5)
mtext("Optimal KrakenUniq kmers = 900, optimal Kraken2 minimizers = 1700")


print(paste0("kmers_max = ",names(rowMeans(F1_matrix_ku))[rowMeans(F1_matrix_ku)==max(rowMeans(F1_matrix_ku))][1],
             ", F1_kmers_max = ",max(rowMeans(F1_matrix_ku))," +/- ",
             sd(F1_matrix_ku[rowMeans(F1_matrix_ku)==max(rowMeans(F1_matrix_ku)),][1,])))
#print(paste0("minimizers_max = ",names(rowMeans(F1_matrix_k2))[rowMeans(F1_matrix_k2)==max(rowMeans(F1_matrix_k2))][1],
#             ", F1_minimizers_max = ",max(rowMeans(F1_matrix_k2))," +/- ",
#             sd(F1_matrix_k2[rowMeans(F1_matrix_k2)==max(rowMeans(F1_matrix_k2)),][1,])))
print(paste0("minimizers_max = ",names(rowMeans(F1_matrix_k2))[rowMeans(F1_matrix_k2)==max(rowMeans(F1_matrix_k2))][1],
             ", F1_minimizers_max = ",max(rowMeans(F1_matrix_k2))," +/- ",
             sd(F1_matrix_k2[rowMeans(F1_matrix_k2)==max(rowMeans(F1_matrix_k2)),])))


library("ggplot2")

df_F1<-data.frame(F1_MAX=c(0.872188787978262,0.854799160062318,0.716029411764706,0.666926507445083),
                  F1_SD_MAX=c(0.131660607968757,0.124587032332865,0.0933239872153626,0.0850384800879344),
                  CLASSIFIER=c("KrakenUniq","Kraken2","KrakenUniq","Kraken2"),
                  SIMULATION=c("Regular","Regular","Pathogen","Pathogen"))
df_F1$CLASSIFIER <- factor(df_F1$CLASSIFIER, levels = c("KrakenUniq","Kraken2"))
df_F1


ggplot(df_F1, aes(x=CLASSIFIER, y=F1_MAX, group=SIMULATION, color=SIMULATION)) + geom_point(size=5) + 
  geom_errorbar(aes(ymin=F1_MAX-F1_SD_MAX, ymax=F1_MAX+F1_SD_MAX), width=0.2, position=position_dodge(0.05)) + 
  theme(text = element_text(size=20))


#sessionInfo()
#R version 4.2.3 (2023-03-15)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.6 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=sv_SE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=sv_SE.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=sv_SE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=sv_SE.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggplot2_3.4.3      viridis_0.6.3      viridisLite_0.4.2  matrixStats_0.63.0 pheatmap_1.0.12   

#loaded via a namespace (and not attached):
#  [1] rstudioapi_0.14    magrittr_2.0.3     tidyselect_1.2.0   munsell_0.5.0      colorspace_2.1-0   R6_2.5.1           rlang_1.1.6       
#[8] fansi_1.0.5        dplyr_1.1.2        tools_4.2.3        grid_4.2.3         gtable_0.3.4       utf8_1.2.4         cli_3.6.5         
#[15] withr_2.5.0        tibble_3.2.1       lifecycle_1.0.4    gridExtra_2.3      farver_2.1.1       RColorBrewer_1.1-3 vctrs_0.6.5       
#[22] glue_1.6.2         labeling_0.4.2     compiler_4.2.3     pillar_1.9.0       generics_0.1.3     scales_1.3.0       pkgconfig_2.0.3   
