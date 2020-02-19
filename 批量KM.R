setwd("E://脊索瘤预后预测模型//结果//家琦分析结果//随机筛选指标结果//手术入路 + 斜坡部位 + degree_of_resection    蛋白：e_cad_cyto + Ki_67_2_n + VEGFA_cyto VEGFAp值更小//PFS")
traindata<-read.table("trainData2.csv",header = T,sep=",")

KM<-function(x){
  mean<-mean(traindata[,x])
  status<-ifelse(traindata[,x]<mean,0,1)
  KMdiff<-survdiff(Surv(PFS,PFS_event==1)~status,data=traindata)
  p.value<-1-pchisq(KMdiff$chisq,length(KMdiff$n)-1)
  KM<-data.frame("characteristcs"=x,"p_value"=p.value)
return(KM)
}

VarNames<-colnames(traindata)[c(27:76)]
Univar<-lapply(c(27:76),KM)
Univar<-ldply(Univar,data.frame)
