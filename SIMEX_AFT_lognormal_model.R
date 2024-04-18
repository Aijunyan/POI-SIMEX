library(tidyverse)
library(survival)
library(survivalAnalysis)
library(purrr)
library(dplyr);library(plyr);#library(simex)
library(simexaft)
library(stringr)
#https://www.utstat.toronto.edu/~brunner/oldclass/312s19/lectures/312s19WeibullRegressionWithR1.pdf
#library(stargazer)
#stargazer(ff, type="text", keep.stat="all")

#bio <- cellcount %>% 
#  rename(CD3pCD8n_S = bio1)

#simex.fit <- simexaft(formula = formula, data = bio_surv, SIMEXvariable = "v.inten",
# repeated = FALSE, repind = list(), err.mat = err.mat, B = 200,
# lambda = seq(0,2,0.1), extrapolation = "quadratic", dist = "weibull")
#simex.est<-summary(simex.fit)$coef

#UNIVARIATE ANALYSIS 
setwd("C:/Users/aijun/Desktop/Aijun_Uvic/Joint Model/R code/Jag_Sim/HGSC")
source("./my.simexAFT.R")

cellcount = read.csv("cell_count.csv")
# Import survival data
outcomes = read.csv("coeur_outcomes.csv")
# Columns of interest:
# Survival from diagnosis (months)
# Death disease specific = 1(dead from disease), 0 (alive or non-disease specific death)

##interested two variable Progression.time.for.stat , CA125.progression.time.in.months
# and get rid of unknown...
#outcomes = outcomes %>% filter(!(Death.Disease.specific %in% c("unknown", "Unknown"))) %>% mutate(Death.Disease.specific = as.numeric(Death.Disease.specific))
##did not do any exclude at this point

colnames(outcomes)[1] <- "pNum"
# Residual disease

outcomes %>% group_by(Residual.disease, Residual.cat.) %>% count()

outcomes <- outcomes %>% mutate(residual_disease = case_when(Residual.cat. == 0 ~ 0,
                                                             Residual.disease %in% c("optimal", "<1cm", "microscopic") ~ 0,
                                                             Residual.disease %in% c("-", "Unknown", NA) ~ NA_real_,
                                                             TRUE ~ 1))

# Stage
outcomes <- outcomes %>% mutate(stage = case_when(Stage %in% c("unknown", "Unknown") ~ as.numeric(NA),
                                                  TRUE ~ as.numeric(Stage)))


table(outcomes$stage)
#Figo
outcomes <- outcomes %>% mutate(figo = case_when(Figo %in% c("I", "Ia", "Ib", "IB", "Ic") ~ 1,
                                                 Figo %in% c("II", "IIa", "IIA", "IIb", "IIc") ~ 2,
                                                 Figo %in% c("III", "IIIa", "IIIb", "IIIB", "IIIc") ~ 3,
                                                 Figo %in% c("iV", "IV", "IIIc-IV") ~ 4,
                                                 Figo %in% c("-", "unknown", "Unknown") ~ as.numeric(NA)))

#outcomes$figo.stage<-ifelse(outcomes$Figo<=2,0,1)
outcomes$figo.stage<-ifelse(outcomes$stage<=2,0,1)

outcomes$survDiag<-outcomes$Survival.from.Diagnosis..months.
outcomes$prog_stat<-outcomes$Progression.time.for.stat
outcomes$CA125_prog<-outcomes$CA125.progression.time.in.months

outcomes$death<-as.numeric(outcomes$All.deaths)

outcomes <- outcomes %>% mutate(grade = case_when(Tumour.Grade %in% c("1","2") ~ 0,
                                                  Tumour.Grade %in% c("3","4") ~ 1,
                                                  Tumour.Grade %in% c("-", "unknown", "Unknown") ~ as.numeric(NA)))


#Age
outcomes <- outcomes %>% mutate(age = as.numeric(Age.at.Diagnosis))
outcomes$agegt60<-ifelse(outcomes$age>=60,1,0)
outcomes$resid<-outcomes$residual_disease


dictionary_old <- c(
  "CD3nCD8nPD1p_S" = "CD3-PD1+ ~ S",
  "CD3pCD8n_S" = "T cell, CD8-FoxP3- ~ S",
  "CD3pCD8nPD1n_S" = "T cell, CD8-PD1- ~ S",
  "CD3pCD8nPD1p_S" = "T cell, CD8-PD1+ ~ S",
  "CD3pCD8pPD1n_S" = "T cell, CD8+PD1- ~ S",
  "CD3pCD8pPD1p_S" = "T cell, CD8+PD1+ ~ S",
  "CD3pCD8nFoxP3p_S" = "T reg, CD8-FoxP3+ ~ S",
  "CD3pCD8pFoxP3p_S" = "T cell, CD8+FoxP3+ ~ S",
  "CD68pPDL1n_S" = "Macrophage, CD68+PDL1- ~ S",
  "CD68pPDL1p_S" = "Macrophage, CD68+PDL1+ ~ S",
  "CD20p_S" = "B cell, CD79a+CD20+ ~ S",
  "CD79ApCD20n_S" = "Plasma cell, CD79a+CD20- ~ S",
  "CD3nCD8nPD1p_E" = "CD3-PD1+ ~ E",
  "CD3pCD8n_E" = "T cell, CD8-FoxP3- ~ E",
  "CD3pCD8nPD1n_E" = "T cell, CD8-PD1- ~ E",
  "CD3pCD8nPD1p_E" = "T cell, CD8-PD1+ ~ E",
  "CD3pCD8pPD1n_E" = "T cell, CD8+PD1- ~ E",
  "CD3pCD8pPD1p_E" = "T cell, CD8+PD1+ ~ E",
  "CD3pCD8nFoxP3p_E" = "T reg, CD8-FoxP3+ ~ E",
  "CD3pCD8pFoxP3p_E" = "T cell, CD8+FoxP3+ ~ E",
  "CD68pPDL1n_E" = "Macrophage, CD68+PDL1- ~ E",
  "CD68pPDL1p_E" = "Macrophage, CD68+PDL1+ ~ E",
  "CD20p_E" = "B cell, CD79a+CD20+~ E",
  "CD79ApCD20n_E" = "Plasma cell, CD79a+CD20- ~ E",
  "avg_str_prop" = "Stroma percentage (per 1%)",
  "cluster_ID" = "Cluster ID",
  "STR_PROP_CAT" = "High stroma",
  "residual_disease" = "Non-optimal debulking status",
  "stage" = "Stage",
  "age" = "Age"
)

dictionary <- c(
  "CD3pCD8n_S" = "T cell, CD8-FoxP3- ~ S",
  "CD3pCD8nPD1n_S" = "T cell, CD8-PD1- ~ S",
  "CD3pCD8nPD1p_S" = "T cell, CD8-PD1+ ~ S",
  "CD3pCD8pPD1n_S" = "T cell, CD8+PD1- ~ S",
  "CD3pCD8pPD1p_S" = "T cell, CD8+PD1+ ~ S",
  "CD3pCD8nFoxP3p_S" = "T reg, CD8-FoxP3+ ~ S",
  "CD3pCD8pFoxP3p_S" = "T cell, CD8+FoxP3+ ~ S",
  "CD68pPDL1n_S" = "Macrophage, CD68+PDL1- ~ S",
  "CD68pPDL1p_S" = "Macrophage, CD68+PDL1+ ~ S",
  "CD20p_S" = "B cell, CD79a+CD20+ ~ S",
  "CD79ApCD20n_S" = "Plasma cell, CD79a+CD20- ~ S",
  "CD3pCD8n_E" = "T cell, CD8-FoxP3- ~ E",
  "CD3pCD8nPD1n_E" = "T cell, CD8-PD1- ~ E",
  "CD3pCD8nPD1p_E" = "T cell, CD8-PD1+ ~ E",
  "CD3pCD8pPD1n_E" = "T cell, CD8+PD1- ~ E",
  "CD3pCD8pPD1p_E" = "T cell, CD8+PD1+ ~ E",
  "CD3pCD8nFoxP3p_E" = "T reg, CD8-FoxP3+ ~ E",
  "CD3pCD8pFoxP3p_E" = "T cell, CD8+FoxP3+ ~ E",
  "CD68pPDL1n_E" = "Macrophage, CD68+PDL1- ~ E",
  "CD68pPDL1p_E" = "Macrophage, CD68+PDL1+ ~ E",
  "CD20p_E" = "B cell, CD79a+CD20+~ E",
  "CD79ApCD20n_E" = "Plasma cell, CD79a+CD20- ~ E",
  "avg_str_prop" = "Stroma percentage (per 1%)",
  "cluster_ID" = "Cluster ID",
  "STR_PROP_CAT" = "High stroma",
  "residual_disease" = "Non-optimal debulking status",
  "stage" = "Stage",
  "age" = "Age"
)


cellcount <- cellcount %>% filter(pNum.BT %in% outcomes$pNum)
cellcount$pNum<-cellcount$pNum.BT
survdata<-subset(outcomes,select=c("pNum","figo.stage","agegt60","grade","resid","death","survDiag"))
                              
TMA<-substr(cellcount$tma_SeC_row_col,1,1)


TMA.core<-cellcount$tma_SeC_row_col
oldname<-names(cellcount)


##randomly select from one sample
## add one dummy to bio
#bio<-data.frame(bio1=cellcount$CD3pCD8nPD1n_E,bio2=cellcount$CD68pPDL1p_E,A=cellcount$area_mm2_E,pNum=cellcount$pNum.BT,TMA=TMA,TMA.core=TMA.core)

##  1) random select one
#bio_one <- cellcount %>%
#  group_by(pNum) %>%
#  sample_n(size=1)

## 2) pick first one
#bio_one <-cellcount %>%
#   group_by(pNum) %>%
#   slice(1)

##2) compute the total cell count

twomarker<-subset(cellcount,select=c("pNum","CD3pCD8nFoxP3p_E","CD3pCD8nFoxP3p_S","area_mm2_E.BT","area_mm2_S.BT"))


#bio_one<-aggregate(. ~ pNum, twomarker, sum)

## 2) pick first one
bio_one <-twomarker %>%
   group_by(pNum) %>%
   slice(1)

bio_surv<-merge(bio_one,survdata,by="pNum")
bio_surv <- subset(bio_surv,bio_surv$survDiag>0)
bio_surv$logsurv<-log(bio_surv$survDiag)

#write.csv(bio_surv,row.names=F,file="./output/AFT/bio_surv.csv")

#write.csv(NA_ct,row.names=T,file="./output/AFT/bio_surv_NA_ct.csv")

suvdata_bio<-subset(bio_surv,select=c("death","agegt60","figo.stage","resid"))

death_age<-ddply(suvdata_bio, c('death','agegt60'), nrow)

death_figo<-ddply(suvdata_bio, c('death','figo.stage'), nrow)
death_resid<-ddply(suvdata_bio, c('death','resid'), nrow)
colnames(death_age)<-c("death","level","count");colnames(death_figo)<-c("death","level","count");
colnames(death_resid)<-c("death","level","count")

surv_bio_summary<-rbind(cbind(death_age,var="agegt60"), cbind(death_figo,var="figo.stage"),cbind(death_resid,var="resid"))

cell_col=c("CD3pCD8nFoxP3p_E","CD3pCD8nFoxP3p_S")
area_2=c("area_mm2_E.BT","area_mm2_S.BT")

totbios=length(cell_col)
naive.est<-c(); simex.est<-c()
stat<-c()
for (i in 1 :totbios){
    #bio_surv$v.inten<-bio_surv$PDL1pCKp_S/bio_surv$area_mm2_S.BT
    ##this is not right here it is only one variable
    #bio_surv$bio<-bio_surv$CD3pCD8n_E/bio_surv$area_mm2_E.BT
  
   if (str_detect(cell_col[i],paste0('_E','$'))) {
    bio_surv$bio<-bio_surv[,cell_col[i]]/bio_surv[,area_2[1]]
    bio_surv$area<-bio_surv$area_mm2_E.BT
   }
   if (str_detect(cell_col[i],paste0('_S','$'))) {
    bio_surv$bio<-bio_surv[,cell_col[i]]/bio_surv[,area_2[2]]
    bio_surv$area<-bio_surv$area_mm2_S.BT
   }
    
    bio_surv$bio.ct<-bio_surv[,cell_col[i]]
    bio_surv_fm<-subset(bio_surv,select=c("survDiag","logsurv","death","agegt60","resid","figo.stage","bio","area","bio.ct"))

    #err.mat<-sqrt(mean(bio_surv_fm$bio))
   
    #formula<-Surv(survDiag,death) ~ bio+agegt60+grade+resid+figo.stage
    stat_<-summary(bio_surv_fm$bio)
    stat2_<-summary(bio_surv_fm$bio.ct)
    stat_<-cbind(para=cell_col[i],t(stat_),t(stat2_))
    
    formula<-Surv(survDiag,death) ~ bio+agegt60+resid+figo.stage
    naive.fit = survreg(formula=formula,data=bio_surv_fm,dist='lognormal',robust=T)
   


naive.est0<-summary(naive.fit)$table
para<-row.names(naive.est0)

naive.est_<-cbind(naive.est0,para,cell=cell_col[i])

simexaft.est<-POI.simexaft(formula=formula,data=bio_surv_fm,SIMEXvariable="bio",areaVariable="area",repind=list(),B=200,lambda=seq(0,2,0.1),extrapolation="quadratic",dist="lognormal")
simex0<-cbind(simexaft.est[1]$coefficients,simexaft.est[2]$se,simexaft.est[3]$scalereg,
              simexaft.est[4]$pvalue)
colnames(simex0)<-c("coef","se","scale","pvalue")
simex.est_<-cbind(para=row.names(simex0),simex0,cell=cell_col[i])

naive.est<- rbind(naive.est,naive.est_)
simex.est<- rbind(simex.est,simex.est_)
stat<-rbind(stat,stat_)

}
write.csv(naive.est,row.names=F,file="./output/AFT/two_biomarker/stat.csv")
write.csv(naive.est,row.names=F,file="./output/AFT/two_biomarker/naive.est.csv")
write.csv(simex.est,row.names=F,file="./output/AFT/two_biomarker/simex.est.csv")
write.csv(surv_bio_summary,row.names=F,file="./output/AFT/two_biomarker/surv_bio_summary.csv")