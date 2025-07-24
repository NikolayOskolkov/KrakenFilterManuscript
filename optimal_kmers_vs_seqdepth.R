# To run the script and reproduce the figures from the manuscript

# Refining filtering criteria for accurate taxonomic classification of ancient metagenomic data
# Nikolay Oskolkov
# bioRxiv 2025.03.31.646431; doi: https://doi.org/10.1101/2025.03.31.646431 

# please clone the repository https://github.com/NikolayOskolkov/KrakenFilterManuscript


#####################################################################################################################################
###################################################### REGULAR MICROBES #############################################################
#####################################################################################################################################


################################################ GROUND TRUTH MAIN AMETA ANALYSIS ####################################################

library("pheatmap")

#gt_path<-"/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/data_and_lists_for_main_aMeta_analysis/"
gt_path<-"path_to_cloned_github_repo/krakenuniq_sim_regular_microbial/"

ground_truth_matrix<-read.delim(paste0(gt_path,"ground_truth_number_of_reads.txt"),row.names=1,header=TRUE,sep="\t")
rownames(ground_truth_matrix)<-gsub("\\.fa","",rownames(ground_truth_matrix))
ground_truth_matrix<-ground_truth_matrix[order(rownames(ground_truth_matrix)),]
colnames(ground_truth_matrix)<-paste0("Sample",seq(1,10,1))
pheatmap(ground_truth_matrix,display_numbers=TRUE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Ground truth: regular microbial species",number_format="%i")
ground_truth_matrix[1:5,]
colSums(ground_truth_matrix)
seq_depth_regular<-colSums(ground_truth_matrix)

ground_truth_matrix_binary<-ground_truth_matrix
ground_truth_matrix_binary[ground_truth_matrix_binary>0]<-1
ground_truth_matrix_binary[ground_truth_matrix_binary<=0]<-0
pheatmap(ground_truth_matrix_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,legend_breaks = c(0,1),
         main="Binary ground truth: regular microbial species")
ground_truth_matrix_binary[1:5,1:5]
colSums(ground_truth_matrix_binary)


################################################# KRAKENUNIQ MAIN AMETA ANALYSIS ####################################################

#ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
ku_path<-"path_to_cloned_github_repo/krakenuniq_sim_regular_microbial/"

library("matrixStats")
max_kmers<-10000; step_kmers<-10
max_taxReads<-300; step_taxReads<-10
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
taxReads_vector<-seq(from=0,to=max_taxReads,by=step_taxReads)
sample_vector<-seq(from=1,to=10,by=1)
IoU_array<-array(rep(NA, length(sample_vector)*length(taxReads_vector)*length(kmers_vector)),
                 c(length(sample_vector), length(taxReads_vector), length(kmers_vector)))
optimal_kmers_regular<-vector()
optimal_IoU_regular<-vector()
for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  IoU_matrix<-matrix(NA,nrow=length(kmers_vector),ncol=length(taxReads_vector))
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
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
      length(intersect(true_list,query_list))
      IoU_matrix[i,j]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
      
      IoU_array[s,j,i]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
    }
  }
  colnames(IoU_matrix)<-taxReads_vector
  rownames(IoU_matrix)<-kmers_vector
  print(head(IoU_matrix))
  
  filled.contour(IoU_matrix,plot.axes= {
    axis(2,at=as.numeric(colnames(IoU_matrix))/max_taxReads,labels=as.numeric(colnames(IoU_matrix)))
    axis(1,at=as.numeric(rownames(IoU_matrix))/max_kmers,las=2,labels=as.numeric(rownames(IoU_matrix)))},
    nlevels=20,color.palette=terrain.colors,main=paste0("Sample",sample_vector[s]))
  
  mtext(paste0("Optimal kmers = ",rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1],", IoU_max = ",max(IoU_matrix)))
  
  optimal_kmers_regular<-append(optimal_kmers_regular,as.numeric(rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1]))
  optimal_IoU_regular<-append(optimal_IoU_regular,max(IoU_matrix))
  
  print(which(IoU_matrix==max(IoU_matrix), arr.ind=T))
}

optimal_kmers_regular
optimal_IoU_regular
seq_depth_regular

plot(optimal_kmers_regular~seq_depth_regular)


#####################################################################################################################################
###################################################### PATHOGENIC MICROBES ##########################################################
#####################################################################################################################################



################################################ GROUND TRUTH MAIN AMETA ANALYSIS ####################################################

library("pheatmap")

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
seq_depth_pathogenic<-colSums(ground_truth_matrix)

ground_truth_matrix_binary<-ground_truth_matrix
ground_truth_matrix_binary[ground_truth_matrix_binary>0]<-1
ground_truth_matrix_binary[ground_truth_matrix_binary<=0]<-0
pheatmap(ground_truth_matrix_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,legend_breaks = c(0,1),
         main="Binary ground truth: pathogenic microbial species")
ground_truth_matrix_binary[1:5,1:5]
colSums(ground_truth_matrix_binary)

################################################# KRAKENUNIQ MAIN AMETA ANALYSIS ####################################################

#ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_pathogens_unfiltered/"
ku_path<-"path_to_cloned_github_repo/krakenuniq_sim_pathogenic_microbial/"

library("matrixStats")
max_kmers<-10000; step_kmers<-10
max_taxReads<-300; step_taxReads<-10
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
taxReads_vector<-seq(from=0,to=max_taxReads,by=step_taxReads)
sample_vector<-seq(from=1,to=10,by=1)
IoU_array<-array(rep(NA, length(sample_vector)*length(taxReads_vector)*length(kmers_vector)),
                 c(length(sample_vector), length(taxReads_vector), length(kmers_vector)))
optimal_kmers_pathogenic<-vector()
optimal_IoU_pathogenic<-vector()
for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  IoU_matrix<-matrix(NA,nrow=length(kmers_vector),ncol=length(taxReads_vector))
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
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
      length(intersect(true_list,query_list))
      IoU_matrix[i,j]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
      
      IoU_array[s,j,i]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
    }
  }
  colnames(IoU_matrix)<-taxReads_vector
  rownames(IoU_matrix)<-kmers_vector
  print(head(IoU_matrix))
  
  filled.contour(IoU_matrix,plot.axes= {
    axis(2,at=as.numeric(colnames(IoU_matrix))/max_taxReads,labels=as.numeric(colnames(IoU_matrix)))
    axis(1,at=as.numeric(rownames(IoU_matrix))/max_kmers,las=2,labels=as.numeric(rownames(IoU_matrix)))},
    nlevels=20,color.palette=terrain.colors,main=paste0("Sample",sample_vector[s]))
  
  mtext(paste0("Optimal kmers = ",rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1],", IoU_max = ",max(IoU_matrix)))
  
  optimal_kmers_pathogenic<-append(optimal_kmers_pathogenic,as.numeric(rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1]))
  optimal_IoU_pathogenic<-append(optimal_IoU_pathogenic,max(IoU_matrix))
  
  print(which(IoU_matrix==max(IoU_matrix), arr.ind=T))
}

optimal_kmers_pathogenic
optimal_IoU_pathogenic
seq_depth_pathogenic

plot(optimal_kmers_pathogenic~seq_depth_pathogenic)


#####################################################################################################################################
############################################################## SEDADNA ##############################################################
#####################################################################################################################################


################################################ GROUND TRUTH MAIN AMETA ANALYSIS ####################################################

library("pheatmap")

#gt_path<-"/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/data_and_lists_for_16_species_simulation/"
gt_path<-"path_to_cloned_github_repo/krakenuniq_sim_eaDNA/"


ground_truth_matrix<-read.delim(paste0(gt_path,"ground_truth_number_of_reads.txt"),row.names=1,header=TRUE,sep="\t")
ground_truth_matrix<-ground_truth_matrix[order(rownames(ground_truth_matrix)),]
pheatmap(ground_truth_matrix,display_numbers=TRUE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Ground truth: sedaDNA samples",number_format="%i")
ground_truth_matrix[1:5,]
colSums(ground_truth_matrix)
seq_depth_sedadna<-colSums(ground_truth_matrix)

ground_truth_matrix_binary<-ground_truth_matrix
ground_truth_matrix_binary[ground_truth_matrix_binary>0]<-1
ground_truth_matrix_binary[ground_truth_matrix_binary<=0]<-0
pheatmap(ground_truth_matrix_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,legend_breaks = c(0,1),
         main="Binary ground truth: sedaDNA samples")
ground_truth_matrix_binary
colSums(ground_truth_matrix_binary)


################################################# KRAKENUNIQ MAIN AMETA ANALYSIS ####################################################

#ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_sedaDNA_simulations_unfiltered/"
ku_path<-"path_to_cloned_github_repo/krakenuniq_sim_eaDNA/"

library("matrixStats")
max_kmers<-10000; step_kmers<-10
max_taxReads<-300; step_taxReads<-10
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
taxReads_vector<-seq(from=0,to=max_taxReads,by=step_taxReads)
sample_vector<-seq(from=1,to=2,by=1)
IoU_array<-array(rep(NA, length(sample_vector)*length(taxReads_vector)*length(kmers_vector)),
                 c(length(sample_vector), length(taxReads_vector), length(kmers_vector)))
optimal_kmers_sedadna<-vector()
optimal_IoU_sedadna<-vector()
for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  IoU_matrix<-matrix(NA,nrow=length(kmers_vector),ncol=length(taxReads_vector))
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
      length(intersect(true_list,query_list))
      IoU_matrix[i,j]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
      
      IoU_array[s,j,i]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
    }
  }
  colnames(IoU_matrix)<-taxReads_vector
  rownames(IoU_matrix)<-kmers_vector
  print(head(IoU_matrix))
  
  filled.contour(IoU_matrix,plot.axes= {
    axis(2,at=as.numeric(colnames(IoU_matrix))/max_taxReads,labels=as.numeric(colnames(IoU_matrix)))
    axis(1,at=as.numeric(rownames(IoU_matrix))/max_kmers,las=2,labels=as.numeric(rownames(IoU_matrix)))},
    nlevels=20,color.palette=terrain.colors,main=paste0("Sample",sample_vector[s]))
  
  mtext(paste0("Optimal kmers = ",rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1],", IoU_max = ",max(IoU_matrix)))
  
  optimal_kmers_sedadna<-append(optimal_kmers_sedadna,as.numeric(rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1]))
  optimal_IoU_sedadna<-append(optimal_IoU_sedadna,max(IoU_matrix))
  
  print(which(IoU_matrix==max(IoU_matrix), arr.ind=T))
}

optimal_kmers_sedadna
optimal_IoU_sedadna
seq_depth_sedadna

plot(optimal_kmers_sedadna~seq_depth_sedadna)




###########################################################################################################################

plot(optimal_kmers_regular~seq_depth_regular, pch=19, col="blue", cex = 1.5, xlim = c(1e+5, 7e+5), ylim = c(0, 2000),
     xlab="Sequencing depth", ylab="Optimal number of unique kmers")
seq_depth_pathogenic<-seq_depth_pathogenic[-which(optimal_kmers_pathogenic==5800)]
optimal_kmers_pathogenic<-optimal_kmers_pathogenic[-which(optimal_kmers_pathogenic==5800)]
points(optimal_kmers_pathogenic~seq_depth_pathogenic, pch=19, col="darkgreen", cex = 1.5)
points(optimal_kmers_sedadna~seq_depth_sedadna, pch=19, col="orange", cex = 1.5)
legend("topleft", inset=.02, c("Regular microbial dataset","Pathogenic microbial dataset","env/sedaDNA dataset"), 
       fill=c("blue","darkgreen","orange"), cex=1.5)

summary(lm(c(optimal_kmers_regular,optimal_kmers_pathogenic,optimal_kmers_sedadna)~
             0+c(seq_depth_regular,seq_depth_pathogenic,seq_depth_sedadna)))
abline(lm(c(optimal_kmers_regular,optimal_kmers_pathogenic,optimal_kmers_sedadna)~
             0+c(seq_depth_regular,seq_depth_pathogenic,seq_depth_sedadna)), col = "red", lwd = 2)

mtext("p-value of linear regression association = 7.3e-7")





#plot(log10(optimal_kmers_regular+1)~log10(seq_depth_regular+1), pch=19, col="blue", cex = 1.5, xlim = c(5, 5.85), ylim=c(2.5, 4))
#points(log10(optimal_kmers_pathogenic+1)~log10(seq_depth_pathogenic+1), pch=19, col="darkgreen", cex = 1.5)
#points(log10(optimal_kmers_sedadna+1)~log10(seq_depth_sedadna+1), pch=19, col="orange", cex = 1.5)

#summary(lm(c(log10(optimal_kmers_regular+1),log10(optimal_kmers_pathogenic+1),log10(optimal_kmers_sedadna+1))~
#             c(log10(seq_depth_regular+1),log10(seq_depth_pathogenic+1),log10(seq_depth_sedadna+1))))
#abline(lm(c(log10(optimal_kmers_regular+1),log10(optimal_kmers_pathogenic+1),log10(optimal_kmers_sedadna+1))~
#            c(log10(seq_depth_regular+1),log10(seq_depth_pathogenic+1),log10(seq_depth_sedadna+1))), col = "red", lwd = 2)




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



