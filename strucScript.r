library(readr)
library(stringr)
library(jsonlite)
library(data.table)
library(magrittr)
library(lazyeval)
library(plyr)
library(dplyr)
getFlame=function(x){
  tryCatch({
  return(str_extract(x,"G[0-9]{1}S[0-9]{1}"))},
  error = function(e) {""})
}

getTumor=function(x){
  tryCatch({
  pattern="肿块[0-9]{1,}\\.{0,1}[0-9]{0,}×[0-9]{1,}\\.{0,1}[0-9]{0,}|结节[0-9]{1,}\\.{0,1}[0-9]{0,}×[0-9]{1,}\\.{0,1}[0-9]{0,}"
  desc=unlist(str_extract_all(x,pattern))
  digits=ifelse(length(desc)==0,NA,as.numeric(unlist(str_extract_all(desc,pattern="[0-9]{1,}\\.{0,1}[0-9]{0,}"))))
  return(list(Diameter=max(digits),Count=ifelse(length(desc)==0,NA,length(desc))))},
error = function(e) {""})
}

getHard=function(x){
  tryCatch({
  pattern="无肝硬化|无明显肝硬化|混合结节型肝硬化|小结节型肝硬化|大结节型肝硬化|未见肝硬化"
  Hard=str_extract(x,pattern)
  HardResult=str_extract(Hard,"无|混合结节型|大结节型|小结节型|结节型|未见")
  HardResult=ifelse(HardResult=="未见","无",HardResult)
  return(HardResult)},
  error = function(e) {""})
}

getOrg=function(x){
  tryCatch({
  pattern="粗梁|细梁|团片|腺管|巢团|透明细胞|梁索|血管网"
  return(str_extract(x,pattern))},
error = function(e) {""})
}

getView=function(x){
  tryCatch({
  pattern="单结节型无包膜|多结节型无包膜|单结节型部分有包膜|多结节型部分有包膜|单结节型有包膜|多结节型有包膜"
  return(str_extract(x,pattern))},
  error = function(e) {""})
}

getStruc=function(x){
  tryCatch({
  pattern="无假小叶结构|未见假小叶结构|见假小叶结构|呈假小叶结构"
  Stru=str_extract(x,pattern)
  StruResult=str_extract(Stru,"无|未见|呈|见")
  if (is.na(StruResult)==T) {StruResult=StruResult} 
  else{
  if (StruResult=="未见" | StruResult=="无") {StruResult="无"}
  if (StruResult=="呈" | StruResult=="见") {StruResult="有"}
  }
  return(StruResult)},
  error = function(e) {""})
}

getBlock=function(x){
  tryCatch({
  pattern="未见血管癌栓|见.{0,4}血管癌栓|有血管癌栓"
  Block=str_extract(x,pattern)
  BlockResult=str_extract(Block,"未见|见|有")
  if (is.na(BlockResult)==T) {BlockResult=BlockResult}
  else {
  if(BlockResult=="未见") {BlockResult="无"}
  if(BlockResult=="见") {BlockResult="有"}
  }
  return(BlockResult)},
  error = function(e) {""})
}

getPos=function(x){
  tryCatch({
  # pattern="肝左叶|肝右叶|肝尾叶|全肝|特殊肝段|肝中叶|部分肝|胆囊床肝组织|胆囊|"
  #return(str_extract(x,pattern))
    data=str_extract(x,"[\\(,（].{1,}?[\\),）]")
    return(gsub(data,pattern = "\\(|\\)|（|）",replacement = ""))},
  error = function(e) {""})
}

getDis=function(x){
  tryCatch({
  pattern="肝细胞癌|小肝细胞癌|肝内胆管癌|海绵状血管瘤|孤立性坏死结节|局灶性结节性增生|腺癌|肝细胞腺瘤"
  return(str_extract(x,pattern))},
  error = function(e) {""})
}

getGrade=function(x){
  tryCatch({
pattern="高度|中度|低度|Ⅰ～Ⅱ|Ⅱ～Ⅲ|Ⅲ～Ⅳ|Ⅱ|Ⅲ|Ⅳ|Ⅰ"
return(str_extract(x,pattern))},
error = function(e) {""})
}

getMVI=function(x){
  tryCatch({
  return(str_extract(x,"M0|M1|M2"))},
  error = function(e) {""})
}

getFlame=function(x){
  tryCatch({
  return(str_extract(x,"G[0-9]{1}S[0-9]{1}"))},
  error = function(e) {""})
}

getBlood=function(x){
  tryCatch({
  x=gsub(x,pattern = "（",replacement="(") %>% gsub(pattern = "）",replacement=")")
  key4=c("AFP","CA199","HBV","HCV")
  keys=unlist(str_split(toupper(str_replace_all(x,pattern="\\([+,-]*\\)|-",replacement = "")),pattern=",|，|。"))
  realKey=intersect(keys,key4)
  if (length(realKey)==0){return(NA)}
  else{
  realValue=str_replace_all(str_extract_all(x,"\\(.{1,3}\\)")[[1]],pattern = "\\(|\\)", replacement = "")
  result=data.frame(t(realValue))
  colnames(result)=realKey
  return(result)}},
  error = function(e) {""})
}

getDye=function(x){
  tryCatch({
  x=gsub(x,pattern = "（",replacement="(") %>% gsub(pattern = "）",replacement=")")
  key2=c("MASSON","AB","VG","网染","弹力","MASSON染色","AB染色","VG染色","PAS")
  keys=unlist(str_split(toupper(str_replace_all(x,pattern="\\([+,-]*\\)|-",replacement = "")),pattern=",|，|。|："))
  realKey=intersect(keys,key2)
  if (length(realKey)==0){return(NA)}
  else{
  realValue=str_replace_all(str_extract_all(x,"\\([+,-]{1,3}\\)")[[1]],pattern = "\\(|\\)", replacement = "")
  result=data.frame(t(realValue))
  colnames(result)=realKey
  return(result)}},
  error = function(e) {""})
  
}

getGene=function(x){
  tryCatch({
  x=gsub(x,pattern = "（",replacement="(") %>% gsub(pattern = "）",replacement=")")
  key3=c("ALK","BCAT","BRAF","CKIT","EGFR","KRAS","HER2","NRAS","PDGFRA","PIK3CA","RET","ROS1","检测号","检测项目")
  pattern=paste0(key3,collapse = "|")
  modified=toupper(str_replace_all(x,pattern="\\([+,-]*\\)|-",replacement = ""))
  values=unlist(str_split(modified,pattern))
  realKey=unlist(str_extract_all(modified,pattern))
  if (length(realKey)==0){return(NA)}
  else{
  result0=data.frame(raw=values[-1])
  result0$value=sapply(result0$raw,function(x) {if(grepl("未",x)){return(FALSE)}; if (grepl("存在",x)){return(TRUE)}})
  result0$key=realKey
  result0=data.table(result0)
  result1=result0[,max(values),by=key]
  result1$mut=ifelse(result1$V1==1,"突变","未突变")
  result=data.frame(t(result1$mut))
  colnames(result)=result1$key
  return(result)}},
  error = function(e) {""})
}

getIHC=function(x){
  tryCatch({
  # key1=c("4EBP1","A103","A1AT","AB","ACE2","AE1/AE3","AFP","AKP","AKR1","AKR1B10","AKT","ALK","ALK\\[克隆号D5F3,VENTANA\\]","ALK1","ALKD5F3","ALKN","Amyloid",
  #        "ANTICRP","APQ1",
  #        "AQP1","AQPI","AR","ARG1","ARGINASE","ARID1","ARID1A","ATRX","BCANTENI","BCAT","BCATENIN","BCATONIN","BCL2","BCL6","BEREP4","BL56","BRCA1","C4D","CA199",
  #        "CA242","CA9","CAIX","CALDESMON","CALPONIN","CAM5.2","CD10","CD105","CD117","CD133","CD138","CD15","CD1A","CD20","CD21","CD23","CD24",
  #        "CD3","CD31","CD34","CD35","CD4","CD44","CD45RO","CD5","CD56","CD61","CD68","CD79<U+2C6D>","CD79A","CD8","CD90","CDK4","CDX2",
  #        "CEA","CERBB2","CGA","CHG","CHRA","CHROA","CHROMOGRANINA","CK","CK14","CK18","CK19","CK20","CK5/6","CK7","CK8","CK广","CMET",
  #        "CMV","CMYC","COX","COX2","CR","CRP","CYCLIN","CYCLIND1","CYCLIND2","D1","D1J2","D240","D2J2","DAXX","DES","DESMIN","DOG1","DPD",
  #        "EBER","ECAD","EGFR","EGFRE746","EGFRL858","EMA","EPCAM","ER","ERCC1","F8","FⅧ","FGF2","FGF218","FGFR1","FIBRINOGEN","FLEX","FLXE",
  #        "FR2","FR34","GAL4","GALECTIN3","GATA3","GFAP","GLY3","GPC3","GRANB","GS","HBCAG","HBME1","HBSAB","HBSAG","HCA\\{BBG066\\}","HCA\\{BCE075\\}",
  #        "HCA\\{腹水\\}","HCK","HCV","HEP1","HEPA","HER2","HMA45","HMB45","HNF4A","HNF4Α","HSP70","HTERT","HVA\\{BBG066\\}","IDH1","IGF2",
  #        "IGG\\(H\\+L\\)",
  #        "INHIBIN","INSULIN","ISL","JVⅡ","JVI","JVII","KAPPA","KI67","KIAA","KIAA0101","KLAA","KP1","LAMBDA","LCA","LFABP","LYSOZYME","M5","M6","MASSON","MAT1","MCA",
  #        "MCH1","MCU1","MDR1","MITF","MLH1","MOC31","MPO","MRP1","MRP14","MRP3","MSH2","MSH6","MTOR","MUC1","MUC2","MUC4","MUC5","MUC5AC","MUC6","MUM1",
  #        "MYO",
  #        "MYOD1","MYOGLOBIN","NAPSINA","NCAD","NDRG1","NES","NF","NKX6","NM23","NSE","NUT","OCT4","OPN","OV6","P16","P27","P40","P4EBP1","P504S",
  #        "P53","P62","P63","P70S","PAKT","PAX8","PCAD","PCEA","PCNA","PD1","PDL1","PDL1(142)","PDL1(288)","PDL1(SP142)","PDL1[142]","PDL1[288]",
  #        "PDL1[BP6]","PDL1[BP6001]","PDL1[E1L]","PDL1[E1L3N]","PDL1[SP142]","PDX1","PER","PERFORIN","PGM1","PKM2","PMH6","PMS1","PMS2","PMTOR",
  #        "PNL2","PR","PS6","PSA","PSAP","PTEN","P糖蛋白","RB","ROS1","S100","S100P","SAA","SHDA","SHDB","SM51","SMA","SMAD4","SOX10","SOX9",
  #        "SPA","SSR2","SSR5","SUOX","SURVIVIN","SUVIVIN","SYN","TFF3","TG","TGFΒ1","TIA","TIA1","TOPIIA","TOPOⅡ","TOPOIIA","TP","TPO","TRIM35",
  #        "TS","TSC2","TTF1","VEGF","VEGFR1","VEGFR2","VEGFR3","VENTANA","VI","VILLIN","VIM","VIMENTIN","WT1","Α1AT","ΑAT","ΒCATENIN",
  #        "ΒTUBULIN","抗酸","六胺银","日期序号","样本号")
  x=gsub(x,pattern = "（",replacement="(") %>% gsub(pattern = "）",replacement=")")
  mat=str_extract_all(x,".{1,}?\\(.{1,}?\\)")[[1]] %>% toupper() %>% str_split(pattern="\\(|\\)",simplify = T)
  if (is.na(mat)==T) {result=NA}
  else{
  realKey=gsub(mat[,1],pattern=":|,|-|：|，",replacement="")
  realValue=mat[,2]
  #keys=str_trim(unlist(str_split(toupper(str_replace_all(x,pattern="\\(.{1,}?\\)|-",replacement = "")),pattern=",|，|。|:|：")))
  #realKey=intersect(keys,key1)
  #realValue=str_replace_all(str_extract_all(x,"\\(.{1,}?\\)")[[1]],pattern = "\\(|\\)", replacement = "")
  result=data.frame(t(realValue))
  colnames(result)=realKey
  }
  return(result)},
  error = function(e) {""})
  
}

########## small funcitons for the single sample ,followings are for actual layout of the system

getEye=function(x){
  tryCatch({
    if (is.na(x)==T) {return(NA)}
    else {
      result=c(getTumor(x)$Diameter,getTumor(x)$Count,getHard(x))
      if (length(unique(as.character(result)))==1) {return(NA)}
      else{
        result=data.frame(t(result))
        colnames(result)=c("肿瘤直径","肿瘤数目","肝硬化")
        return(as.list(result))}}},
    error = function(e) {""})
}

getLen=function(x){
  tryCatch({
    if (is.na(x)==T) {return(NA)}
    else {
      result=c(getOrg(x),getView(x),getStruc(x),getBlock(x))
      if (length(unique(as.character(result)))==1) {return(NA)}
      else{
        result=data.frame(t(result))
        colnames(result)=c("组织类型","肉眼类型","假小叶结构","血管癌栓")
        return(as.list(result))}}},
    error = function(e) {""})
  
}

getConclusion=function(x){
  tryCatch({
    if (is.na(x)==T) {return(NA)}
    else {
      result=c(getPos(x),getDis(x),getGrade(x),getMVI(x),getFlame(x))
      if (length(unique(result))==1) {return(NA)}
      else{
        result=data.frame(t(result))
        colnames(result)=c("患病部位","病理类型","分化程度","MVI分级","肝炎分期")
        return(as.list(result))}}},
    error = function(e) {""})
  
}

getPred=function(x){
  tryCatch({
    if (is.na(x)==T) {return(NA)}
    else {
      keys=c("术式", "单发肿瘤","多发肿瘤","肉眼类型","组织类型","分级","卫星灶","脉管侵犯（巨检／手术所见）","微脉管侵犯（显微镜下所见）",
             "MVI提示风险分级","切除面","小胆管癌栓","肝被膜","胆管侵犯","癌周围肝组织","周围神经侵犯","肝硬化","远处转移",
             "肝炎","肝炎程度","纤维化分期","淋巴结胆囊侵犯","淋巴结","另送膈肌","邻近组织侵犯")
      
      pattern=paste0(keys,collapse = "|")
      # modified=toupper(str_replace_all(x,pattern="\\([+,-]*\\)|-",replacement = ""))
      values=gsub(unlist(str_split(x,pattern)),pattern=":|,", replacement=" ")
      realKey=unlist(str_extract_all(x,pattern))
      result=data.frame(t(values[-1]))
      colnames(result)=realKey
      return(as.list(result))}},
    error = function(e) {""})
}


withSample=function(x,func){
  tryCatch({
    #pattern="[0-9]{2}S[0-9]{5}-[0-9]{3}"
    if (is.na(x)==T){return(NA)}
    else{
      pattern="M[0-9]{4}-[0-9]{4}|另一区域"
      samples=unlist(str_split(x,pattern)) 
      SampleNo=unlist(str_extract_all(x,pattern))
      if(length(SampleNo)==0 & length(samples)==1 & str_detect(samples,"\\(.{1}\\)")==F) {return(NA)}
      if (str_detect(samples[1],"\\(.{1}\\)")==T) {SampleNo=c("未知",SampleNo);samples=c("",samples)}
      else {SampleNo=SampleNo;samples=samples}
      result0=sapply(data.frame(sample=samples[-1])$sample,func,simplify = F)
      result=NULL
      for ( i in 1: length(result0)){
        result=rbindlist(list(result, cbind(c(SampleNo[i]),result0[[i]])),fill = T)
      }
      colnames(result)[1]="样本号"
      return(result)}},
    error = function(e) {""})
}

DFGD <- read_csv("E:/pathology/DFGDRaw.csv")
DATA=as.list(rep(0,nrow(DFGD)))
for (i in 1: length(DATA)) {
  cand=DFGD[i,]
  result=list("base_info"=as.list(select(cand,-肉眼所见,-镜下所见,-特殊检查,-病理诊断,-HBV,-HCV,-CA19.9,-AFP)),
              "emr_info"=list(
                "血清检验"=data.frame("样本号"="未知","HBV"=cand$HBV,"HCV"=cand$HCV,"CA19.9"=cand$CA19.9,"AFP"=cand$AFP),
                "肉眼所见"=getEye(cand$肉眼所见),
                "镜下所见"=getLen(cand$镜下所见),
                "免疫组化"=withSample(cand$特殊检查,getIHC),
                "诊断结论"=getConclusion(cand$病理诊断),
                "特殊染色"=withSample(cand$镜下所见,getDye)))
  DATA[[i]]=result
  setTxtProgressBar(txtProgressBar(min=0,max=1,style = 3),value=i/length(DATA))
  
}


system.time(m=toJSON(DATA,auto_unbox = T))
system.time(write(m,"D:/DFGD1.json"))

###############################

getIHCwithSample=function(x){
  
  tryCatch({
  #pattern="[0-9]{2}S[0-9]{5}-[0-9]{3}"
  if (is.na(x)==T){return(NA)}
  else{
  pattern="M[0-9]{4}-[0-9]{4}|另一区域"
  samples=unlist(str_split(x,pattern)) 
  SampleNo=unlist(str_extract_all(x,pattern))
  if(length(SampleNo)==0 & length(samples)==1 & str_detect(samples,"\\(.{1}\\)")==F) {return(NA)}
  if (str_detect(samples[1],"\\(.{1}\\)")==T) {SampleNo=c("未知",SampleNo);samples=c("",samples)}
  else {SampleNo=SampleNo;samples=samples}
  result0=sapply(data.frame(sample=samples[-1])$sample,getIHC,simplify = F)
  result=NULL
  for ( i in 1: length(result0)){
    result=rbindlist(list(result, cbind(c(SampleNo[i]),result0[[i]])),fill = T)
  }
  colnames(result)[1]="样本号"
  # index=paste("样本",1:nrow(result),sep = "")
  # LIST=list("init"=1)
  # for(i in 1:length(index)){
  #   LIST[[i]]=result[i]
  # }
  # names(LIST)=paste("样本",1:nrow(result),sep = "")
  # return(LIST)
  return(result)}},
  error = function(e) {""})
                  
}

getGenewithSample=function(x){
  if (is.na(x)==T){return(NA)}
  else{
  pattern="[0-9]{2}S[0-9]{5}-[0-9]{3}"
  samples=unlist(str_split(x,pattern))
  SampleNo=unlist(str_extract_all(x,pattern))
  if(length(SampleNo)==0 & length(samples)==1 & str_detect(samples,"\\(.{1}\\)")==F) {return(NA)}
  result0=sapply(data.frame(sample=samples[-1])$sample,getGene,simplify = F)
  result=NULL
  for ( i in 1: length(result0)){
    result=rbindlist(list(result, cbind(c(SampleNo[i]),result0[[i]])),fill = T)
  }
  colnames(result)[1]="样本号"
  # index=paste("样本",1:nrow(result),sep = "")
  # LIST=list("init"=1)
  # for(i in 1:length(index)){
  #   LIST[[i]]=result[i]
  # }
  # names(LIST)=paste("样本",1:nrow(result),sep = "")
  # return(LIST)
  return(result)}
}

getDyewithSample=function(x){
  tryCatch({
  
  if (is.na(x)==T){return(NA)}
  else{
  # pattern="[0-9]{2}S[0-9]{5}-[0-9]{3}"
  pattern="M[0-9]{4}-[0-9]{4}"
  samples=unlist(str_split(x,pattern))
  SampleNo=unlist(str_extract_all(x,pattern))
  if(length(SampleNo)==0 & length(samples)==1 & str_detect(samples,"\\(.{1}\\)")==F) {return(NA)}
  if (str_detect(samples[1],"\\(.{1}\\)")==T) {SampleNo=c("未知",SampleNo);samples=c("",samples)}
  else {SampleNo=SampleNo;samples=samples}
  #SampleNo=ifelse(length(SampleNo0)!=length(samples),c("未知",SampleNo0),SampleNo0)
  result0=sapply(data.frame(sample=samples[-1])$sample,getDye,simplify = F)
  result=NULL
  for ( i in 1: length(result0)){
    result=rbindlist(list(result, cbind(c(SampleNo[i]),result0[[i]])),fill = T)
  }
  colnames(result)[1]="样本号"
  # index=paste("样本",1:nrow(result),sep = "")
  # LIST=list("init"=1)
  # for(i in 1:length(index)){
  #   LIST[[i]]=result[i]
  # }
  # names(LIST)=paste("样本",1:nrow(result),sep = "")
  # return(LIST)
  return(result)}},
  error = function(e) {""})
}

getBloodwithSample=function(x){
  if (is.na(x)==T){return(NA)}
  else{
  pattern="[0-9]{2}S[0-9]{5}-[0-9]{3}"
  samples=unlist(str_split(x,pattern))
  SampleNo0=unlist(str_extract_all(x,pattern))
  SampleNo=ifelse(length(SampleNo0)!=length(samples),c("未知",SampleNo0),SampleNo0)
  result0=sapply(data.frame(sample=samples)$sample,getBlood,simplify = F)
  result=NULL
  for ( i in 1: length(result0)){
    result=rbindlist(list(result, cbind(c(SampleNo[i]),result0[[i]])),fill = T)
  }
  colnames(result)[1]="样本号"
  # index=paste("样本",1:nrow(result),sep = "")
  # LIST=list("init"=1)
  # for(i in 1:length(index)){
  #   LIST[[i]]=result[i]
  # }
  # names(LIST)=paste("样本",1:nrow(result),sep = "")
  # return(LIST)
  return(result)}
}



DFGD <- read_csv("E:/pathology/DFGDRaw.csv")
DATA=as.list(rep(0,nrow(DFGD)))
system.time(for (i in 1: length(DATA)) {
  cand=DFGD[i,]
  result=list("base_info"=as.list(select(cand,-肉眼所见,-镜下所见,-特殊检查,-病理诊断)),
              "emr_info"=list(
       "肉眼所见"=getEye(cand$肉眼所见),
       "镜下所见"=getLen(cand$镜下所见),
       "免疫组化"=getIHCwithSample(cand$特殊检查),
       "病理诊断"=getConclusion(cand$病理诊断),
       "特殊染色"=getDyewithSample(cand$镜下所见)))
  DATA[[i]]=result
  setTxtProgressBar(txtProgressBar(min=0,max=1,style = 3),value=i/length(DATA))
  
})

m=toJSON(DATA,auto_unbox = T)
write(m,"D:/DFGD.json")


