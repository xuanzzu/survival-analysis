library(plyr)
library(pROC)
setwd("F:\\radiomics\\meningioma\\patients data\\R语言 logistic建模")
train<-read.csv(file="train.csv")
test<-read.csv(file = "test.csv")
train<-train[,-1]
test<-test[,-1]
#批量单因素logstic回归
Unilogistic.result<-function(x){
  FML<-as.formula(paste0("质地~",x))
  unilogistic<-glm(FML,family = binomial(link = logit),data = train)
  unilogistic.sum<-summary(unilogistic)
  pvalue<-unilogistic.sum$coefficients[2,4]
  unilogistic.result<-data.frame("characteristics"=x,"pvalue"=pvalue)
  return(unilogistic.result)
}
varnames<-colnames(train)[-1]
Unilogsitic<-lapply(varnames,Unilogistic.result)
Unilogsitic<-ldply(Unilogsitic,data.frame)
#筛选p小于0.05的变量
f<-Unilogsitic$characteristics[Unilogsitic$pvalue<0.05]
#多因素logistic回归
FML<-as.formula(paste0("质地~",paste0(f,collapse = "+")))
Multilogistic<-glm(FML,family=binomial(link = logit),data=train)
logit.step<-step(Multilogistic,direction = "backward")
summary(logit.step)
#最终logistic模型
logit.final<-glm(质地 ~ T1c_wavelet.HHH_firstorder_Maximum + T2_log.sigma.3.0.mm.3D_firstorder_Mean + 
                     T2_wavelet.LHL_firstorder_Median + T2_wavelet.HLL_glcm_Correlation + 
                     T2flair_original_glrlm_LongRunHighGrayLevelEmphasis + ADC_wavelet.LHL_firstorder_Skewness + 
                     ADC_wavelet.HLL_firstorder_Skewness + ADC_wavelet.HLH_firstorder_Median + 
                     ADC_wavelet.HLH_firstorder_Skewness, family = binomial(link = logit), 
                   data = train)
#AUC计算
prob<-predict(logit.final,type = c("response"))
train$prob<-prob
AUC.train<-auc(train$质地,train$prob)

prob.test<-predict(logit.final,type = c("response"),newdata = test)
test$prob<-prob.test
AUC.test<-auc(test$质地,test$prob)
#ROC曲线绘制
ROC<-roc(质地~prob,data = train)
plot(ROC)
