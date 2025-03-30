################################################## F1 SCORE ####################################################

library("ggplot2")

df_F1<-data.frame(F1_MAX=c(0.8721888,0.8523335,0.7160294,0.6433779,0.8001231,0.6007739,0.7055882,0.385375,0.8473241,
                           0.6036845,0.6817794,0.3998581,0.8479263,0.6727273,0.8479263,0.2709091,0.6190476,0.2745826),
                  F1_SD_MAX=c(0.1316606,0.1283951,0.09332399,0.1182888,0.13093,0.1141565,0.1019147,0.07667044,0.1080268,
                              0.09539046,0.1207271,0.08885505,0.01303423,0.1799908,0.01303423,0.1259936,0.0673435,0.1311886),
                  FILTER=c("K","C","K","C","R","K / R","R","K / R","(K / R) * C","(K / R) * dexp(C)",
                           "(K / R) * C","(K / R) * dexp(C)","K","C","R","K / R","(K / R) * C","(K / R) * dexp(C)"),
                  SIMULATION=c("Regular","Regular","Pathogen","Pathogen","Regular","Regular","Pathogen","Pathogen","Regular",
                               "Regular","Pathogen","Pathogen","env/sedaDNA","env/sedaDNA","env/sedaDNA","env/sedaDNA",
                               "env/sedaDNA","env/sedaDNA"))
df_F1$FILTER <- factor(df_F1$FILTER, levels = c("K","C","R","K / R","(K / R) * C","(K / R) * dexp(C)"))
df_F1$SIMULATION <- factor(df_F1$SIMULATION, levels = c("Regular","Pathogen","env/sedaDNA"))
df_F1


ggplot(df_F1, aes(x=FILTER, y=F1_MAX, group=SIMULATION, color=SIMULATION)) + geom_point(size=5) + 
  geom_errorbar(aes(ymin=F1_MAX-F1_SD_MAX, ymax=F1_MAX+F1_SD_MAX), width=0.2, position=position_dodge(0.05)) + 
  theme(text = element_text(size=20))


#########################################################################################################################

df_F1_no_group<-df_F1
df_F1_no_group$SIMULATION<-NULL
df_F1_no_group

subset1<-df_F1_no_group[as.character(df_F1_no_group$FILTER)=="K",]
subset1_mod<-data.frame(F1_MAX=mean(subset1$F1_MAX), F1_SD_MAX=sqrt(sum(subset1$F1_SD_MAX^2)), FILTER="K")
subset1_mod

subset2<-df_F1_no_group[as.character(df_F1_no_group$FILTER)=="C",]
subset2_mod<-data.frame(F1_MAX=mean(subset2$F1_MAX), F1_SD_MAX=sqrt(sum(subset2$F1_SD_MAX^2)), FILTER="C")
subset2_mod

subset3<-df_F1_no_group[as.character(df_F1_no_group$FILTER)=="R",]
subset3_mod<-data.frame(F1_MAX=mean(subset3$F1_MAX), F1_SD_MAX=sqrt(sum(subset3$F1_SD_MAX^2)), FILTER="R")
subset3_mod

subset4<-df_F1_no_group[as.character(df_F1_no_group$FILTER)=="K / R",]
subset4_mod<-data.frame(F1_MAX=mean(subset4$F1_MAX), F1_SD_MAX=sqrt(sum(subset4$F1_SD_MAX^2)), FILTER="K / R")
subset4_mod

subset5<-df_F1_no_group[as.character(df_F1_no_group$FILTER)=="(K / R) * C",]
subset5_mod<-data.frame(F1_MAX=mean(subset5$F1_MAX), F1_SD_MAX=sqrt(sum(subset5$F1_SD_MAX^2)), FILTER="(K / R) * C")
subset5_mod

subset6<-df_F1_no_group[as.character(df_F1_no_group$FILTER)=="(K / R) * dexp(C)",]
subset6_mod<-data.frame(F1_MAX=mean(subset6$F1_MAX), F1_SD_MAX=sqrt(sum(subset6$F1_SD_MAX^2)), FILTER="(K / R) * dexp(C)")
subset6_mod


df_F1<-rbind(rbind(rbind(rbind(rbind(subset1_mod,subset2_mod),subset3_mod),subset4_mod),subset5_mod),subset6_mod)
df_F1$FILTER <- factor(df_F1$FILTER, levels = c("K","C","R","K / R","(K / R) * C","(K / R) * dexp(C)"))
df_F1

ggplot(df_F1, aes(x=FILTER, y=F1_MAX)) + geom_point(size=5) + 
  geom_errorbar(aes(ymin=F1_MAX-F1_SD_MAX, ymax=F1_MAX+F1_SD_MAX), width=0.2, 
                position=position_dodge(0.05)) + theme(text = element_text(size=20))













############################################ JACCARD SIMILARITY ################################################

library("ggplot2")

df_JACCARD<-data.frame(JACCARD_MAX=c(0.5645455,0.484697,0.7910196,0.7585985,0.6823123,0.4379891,0.5531935,0.2412924,0.746549,
                                     0.4381647,0.5289394,0.2534114,0.73611111,0.52083333,0.7361111,0.1597561,0.45,0.1625),
                  JACCARD_SD_MAX=c(0.1054581,0.1332891,0.1689929,0.1586859,0.1590277,0.1177216,0.1146854,0.06127184,0.1361275,
                                   0.09488602,0.1430131,0.07056223,0.01964186,0.2062395,0.01964186,0.08450788,0.07071068,
                                   0.08838835),
                  FILTER=c("K","C","K","C","R","K / R","R","K / R","(K / R) * C","(K / R) * dexp(C)",
                           "(K / R) * C","(K / R) * dexp(C)","K","C","R","K / R","(K / R) * C","(K / R) * dexp(C)"),
                  SIMULATION=c("Pathogen","Pathogen","Regular","Regular","Regular","Regular","Pathogen","Pathogen","Regular",
                               "Regular","Pathogen","Pathogen",
                               "env/sedaDNA","env/sedaDNA","env/sedaDNA","env/sedaDNA","env/sedaDNA","env/sedaDNA"))
df_JACCARD$FILTER <- factor(df_JACCARD$FILTER, levels = c("K","C","R","K / R","(K / R) * C","(K / R) * dexp(C)"))
df_JACCARD$SIMULATION <- factor(df_JACCARD$SIMULATION, levels = c("Regular","Pathogen","env/sedaDNA"))
df_JACCARD


ggplot(df_JACCARD, aes(x=FILTER, y=JACCARD_MAX, group=SIMULATION, color=SIMULATION)) + geom_point(size=5) + 
  geom_errorbar(aes(ymin=JACCARD_MAX-JACCARD_SD_MAX, ymax=JACCARD_MAX+JACCARD_SD_MAX), width=0.2, 
                position=position_dodge(0.05)) + theme(text = element_text(size=20))


#########################################################################################################################

df_JACCARD_no_group<-df_JACCARD
df_JACCARD_no_group$SIMULATION<-NULL
df_JACCARD_no_group

subset1<-df_JACCARD_no_group[as.character(df_JACCARD_no_group$FILTER)=="K",]
subset1_mod<-data.frame(JACCARD_MAX=mean(subset1$JACCARD_MAX), JACCARD_SD_MAX=sqrt(sum(subset1$JACCARD_SD_MAX^2)), FILTER="K")
subset1_mod

subset2<-df_JACCARD_no_group[as.character(df_JACCARD_no_group$FILTER)=="C",]
subset2_mod<-data.frame(JACCARD_MAX=mean(subset2$JACCARD_MAX), JACCARD_SD_MAX=sqrt(sum(subset2$JACCARD_SD_MAX^2)), FILTER="C")
subset2_mod

subset3<-df_JACCARD_no_group[as.character(df_JACCARD_no_group$FILTER)=="R",]
subset3_mod<-data.frame(JACCARD_MAX=mean(subset3$JACCARD_MAX), JACCARD_SD_MAX=sqrt(sum(subset3$JACCARD_SD_MAX^2)), FILTER="R")
subset3_mod

subset4<-df_JACCARD_no_group[as.character(df_JACCARD_no_group$FILTER)=="K / R",]
subset4_mod<-data.frame(JACCARD_MAX=mean(subset4$JACCARD_MAX), JACCARD_SD_MAX=sqrt(sum(subset4$JACCARD_SD_MAX^2)), FILTER="K / R")
subset4_mod

subset5<-df_JACCARD_no_group[as.character(df_JACCARD_no_group$FILTER)=="(K / R) * C",]
subset5_mod<-data.frame(JACCARD_MAX=mean(subset5$JACCARD_MAX), JACCARD_SD_MAX=sqrt(sum(subset5$JACCARD_SD_MAX^2)), FILTER="(K / R) * C")
subset5_mod

subset6<-df_JACCARD_no_group[as.character(df_JACCARD_no_group$FILTER)=="(K / R) * dexp(C)",]
subset6_mod<-data.frame(JACCARD_MAX=mean(subset6$JACCARD_MAX), JACCARD_SD_MAX=sqrt(sum(subset6$JACCARD_SD_MAX^2)), FILTER="(K / R) * dexp(C)")
subset6_mod


df_JACCARD<-rbind(rbind(rbind(rbind(rbind(subset1_mod,subset2_mod),subset3_mod),subset4_mod),subset5_mod),subset6_mod)
df_JACCARD$FILTER <- factor(df_JACCARD$FILTER, levels = c("K","C","R","K / R","(K / R) * C","(K / R) * dexp(C)"))
df_JACCARD

ggplot(df_JACCARD, aes(x=FILTER, y=JACCARD_MAX)) + geom_point(size=5) + 
  geom_errorbar(aes(ymin=JACCARD_MAX-JACCARD_SD_MAX, ymax=JACCARD_MAX+JACCARD_SD_MAX), width=0.2, 
                position=position_dodge(0.05)) + theme(text = element_text(size=20))
