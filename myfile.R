# myfile.R


#* @get /IHC

singleIHC=function(input){
library(readxl)
library(plyr)
library(dplyr)
library(stringr)
library(data.table)
library(lazyeval)
library(jsonlite)
# key1=c("4EBP1","A103","A1AT","AB","ACE2","AE1/AE3","AFP","AKP","AKR1","AKR1B10","AKT","ALK","ALK[克隆号D5F3,VENTANA]","ALK1","ALKD5F3","ALKN","Amyloid",
#        "ANTICRP","APQ1","AQP1","AQPI","AR","ARG1","ARGINASE","ARID1","ARID1A","BCANTENI","BCAT","BCATENIN","BCATONIN","BCL2","BCL6","BEREP4","BL56","BRCA1","C4D","CA199",
#        "CA242","CA9","CAIX","CALDESMON","CALPONIN","CAM5.2","CD10","CD105","CD117","CD133","CD138","CD15","CD1A","CD20","CD21","CD23","CD24","CD3","CD31","CD34","CD35","CD4","CD44","CD45RO","CD5","CD56","CD61","CD68","CD79<U+2C6D>","CD79A","CD8","CD90","CDK4","CDX2",
#        "CEA","CERBB2","CGA","CHG","CHRA","CHROA","CHROMOGRANINA","CK","CK14","CK18","CK19","CK20","CK5/6","CK7","CK8","CK广","CMET","CMV","CMYC","COX","COX2","CR","CRP","CYCLIN","CYCLIND1","CYCLIND2","D1","D1J2","D240","D2J2","DAXX","DES","DESMIN","DOG1","DPD",
#        "EBER","ECAD","EGFR","EGFRE746","EGFRL858","EMA","EPCAM","ER","ERCC1","F8","FⅧ","FGF2","FGF218","FGFR1","FIBRINOGEN","FLEX","FLXE","FR2","FR34","GAL4","GALECTIN3","GATA3","GFAP","GLY3","GPC3","GRANB","GS","HBCAG","HBME1","HBSAB","HBSAG","HCA\\{BBG066\\}","HCA\\{BCE075\\}",
#        "HCA\\{腹水\\}","HCK","HCV","HEP1","HEPA","HER2","HMA45","HMB45","HNF4A","HNF4Α","HSP70","HTERT","HVA\\{BBG066\\}","IDH1","IGG\\(H\\+L\\)",
#        "INHIBIN","INSULIN","ISL","JVⅡ","JVI","JVII","KAPPA","KI67","KIAA","KIAA0101","KLAA","KP1","LAMBDA","LCA","LFABP","LYSOZYME","M5","M6","MASSON","MAT1","MCA",
#        "MCH1","MCU1","MDR1","MITF","MLH1","MOC31","MPO","MRP1","MRP14","MRP3","MSH2","MSH6","MTOR","MUC1","MUC2","MUC4","MUC5","MUC5AC","MUC6","MUM1",
#        "MYOD1","MYOGLOBIN","NAPSINA","NCAD","NDRG1","NES","NF","NKX6","NM23","NSE","NUT","OCT4","OPN","OV6","P16","P27","P40","P4EBP1","P504S",
#        "P53","P62","P63","P70S","PAKT","PAX8","PCAD","PCEA","PCNA","PD1","PDL1","PDL1(142)","PDL1(288)","PDL1(SP142)","PDL1[142]","PDL1[288]",
#        "PDL1[BP6]","PDL1[BP6001]","PDL1[E1L]","PDL1[E1L3N]","PDL1[SP142]","PDX1","PER","PERFORIN","PGM1","PKM2","PMH6","PMS1","PMS2","PMTOR","PNL2","PR","PS6","PSA","PSAP","PTEN","P糖蛋白","RB","ROS1","S100","S100P","SAA","SHDA","SHDB","SM51","SMA","SMAD4","SOX10","SOX9", "SPA","SSR2","SSR5","SUOX","SURVIVIN","SUVIVIN","SYN","TFF3","TG","TGFΒ1","TIA","TIA1","TOPIIA","TOPOⅡ","TOPOIIA","TP","TPO","TRIM35",
#        "TS","TSC2","TTF1","VEGF","VEGFR1","VEGFR2","VEGFR3","VENTANA","VI","VILLIN","VIM","VIMENTIN","WT1","Α1AT","ΑAT","ΒCATENIN",
#        "ΒTUBULIN","抗酸","六胺银","日期序号","样本号")
key1=as.character(read.table("E:/pathology/dictOrdered.txt")$V1)
getIHC=function(x,y){
  keys=unlist(str_split(toupper(str_replace_all(x,pattern="\\([+,-]*\\)|-",replacement = "")),pattern=","))
  index=which(keys==toupper(y))
  real=unlist(str_split(x,pattern=","))[index]
  return(unlist(str_split(real,pattern = "\\(|\\)"))[2])
}

IHC=NULL

for (i in 1: length(key1)){

  tryCatch (
    {
      result=as.character(sapply(input,getIHC,y=key1[i]))
      if (result=="NULL"){
        IHC=IHC
      }
      else {
        IHC=cbind(IHC,c(key1[i],result))
      }

    },
    error = function(e) {""})
}
  colnames(IHC)=IHC[1,]

  return(toJSON(data.frame(IHC)[-1,]))
}

