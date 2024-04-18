library(tidyverse)
library(survival)
library(survivalAnalysis)
library(purrr)
library(dplyr);library(plyr);#library(simex)
library(simexaft)
library(stringr)
library(boot);library(sjstats)

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
bio_one <-cellcount %>%
   group_by(pNum) %>%
   slice(1)

bio_surv<-merge(bio_one,survdata,by="pNum")

NA_ct<-sapply(bio_surv, function(x) sum(is.na(x)))

write.csv(bio_surv,row.names=F,file="./output/AFT/bio_surv.csv")

write.csv(NA_ct,row.names=T,file="./output/AFT/bio_surv_NA_ct.csv")

#CD3pCD8nPD1n_E, CD68pPDL1p_E


#$ endwith ^startwith
cell_S<-names(bio_one)[str_detect(names(bio_one),paste0('_S','$'))]
cell_E<-names(bio_one)[str_detect(names(bio_one),paste0('_E','$'))]

#[1] "CD3pCD8n_S"       "CD3pCD8nFoxP3p_S" "CD3pCD8p_S"       "CD20p_S"         
#[5] "CD3pCD8pFoxP3p_S" "CD79ApCD20n_S"    "CD3pCD8nPD1n_S"   "CD3pCD8nPD1p_S"  
#[9] "CD3pCD8pPD1n_S"   "CD3pCD8pPD1p_S"   "CD68nPDL1nPD1p_S" "area_mm2_S"      
#[13] "CD68pPDL1n_S"     "CD68pPDL1p_S"     "PDL1pCKp_S"      
#[1] "CD3pCD8n_E"       "CD3pCD8nFoxP3p_E" "CD3pCD8p_E"       "CD20p_E"         
#[5] "CD3pCD8pFoxP3p_E" "CD79ApCD20n_E"    "CD3pCD8nPD1n_E"   "CD3pCD8nPD1p_E"  
#[9] "CD3pCD8pPD1n_E"   "CD3pCD8pPD1p_E"   "CD68nPDL1nPD1p_E" "area_mm2_E"      
#[13] "CD68pPDL1n_E"     "CD68pPDL1p_E"     "PDL1pCKp_E"

##START 8  then a lot NA
#focus on those without NA cell cout;

cell_col=c("CD3pCD8n_E","CD3pCD8n_S","CD3pCD8nFoxP3p_E","CD3pCD8nFoxP3p_S","CD3pCD8p_E","CD3pCD8p_S",
           "CD20p_E","CD20p_S","CD3pCD8pFoxP3p_E","CD3pCD8pFoxP3p_S","CD79ApCD20n_E","CD79ApCD20n_S")
area_2=c("area_mm2_E.BT","area_mm2_S.BT")

totbios=length(cell_col)
naive.est<-c(); simex.est<-c()

    #bio_surv$v.inten<-bio_surv$PDL1pCKp_S/bio_surv$area_mm2_S.BT
    ##this is not right here it is only one variable
    #bio_surv$bio<-bio_surv$CD3pCD8n_E/bio_surv$area_mm2_E.BT

   i<-1
   if (str_detect(cell_col[i],paste0('_E','$'))) {
    bio_surv$bio<-bio_surv[,cell_col[i]]/bio_surv[,area_2[1]]
    bio_surv$area<-bio_surv$area_mm2_E.BT
   }
   
   if (str_detect(cell_col[i],paste0('_S','$')))   {
    bio_surv$bio<-bio_surv[,cell_col[i]]/bio_surv[,area_2[2]]
    bio_surv$area<-bio_surv$area_mm2_S.BT
   }
  
    bio_surv_fm<-subset(bio_surv,select=c("survDiag","death","agegt60","resid","figo.stage","bio","area"))
  
aft.coef<-function(data,ind){  
    #err.mat<-sqrt(mean(bio_surv_fm$bio))
    dat<-data[ind,]
    #formula<-Surv(survDiag,death) ~ bio+agegt60+grade+resid+figo.stage
    formula<-Surv(survDiag,death) ~ bio+agegt60+resid+figo.stage
   simexaft.est<-POI.simexaft(formula=formula,data=dat,SIMEXvariable="bio",areaVariable="area",repind=list(),B=200,lambda=seq(0,2,0.1),extrapolation="quadratic",dist="weibull")
   simex0<-t(simexaft.est[1]$coefficients)
   
   simex.est<-cbind(para=row.names(simex0),simex0,cell=cell_col[i])
   
   #simex.est<-subset(simex0,select="bio")
return (simex.est)
}

test<-sample(bio_surv_fm, nrow(bio_surv_fm), replace = T, prob = NULL)



