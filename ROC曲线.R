setwd("E://脊索瘤预后预测模型//结果//家琦分析结果//随机筛选指标结果//手术入路 + 斜坡部位 + degree_of_resection    蛋白：e_cad_cyto + Ki_67_2_n + VEGFA_cyto VEGFAp值更小//PFS")
library(survival)
library(timeROC)
traindata<-read.table("trainData.csv",header = T,sep=",")
testdata<-read.table("testData.csv",header = T,sep=",")
y<-Surv(traindata$PFS,traindata$PFS_event==1)
MultiCox.final<-coxph(y~ degree_of_resection + e_cad_cyto + Ki_67_2_n + VEGFA_cyto,ties = c("breslow"),
                      data = traindata)
traindata$lp<-predict(MultiCox.final,type = "lp")
testdata$lp<-predict(MultiCox.final,type = "lp",newdata = testdata)
ROC<-timeROC(T=traindata$PFS,delta = traindata$PFS_event,
             marker = traindata$lp,cause = 1,weighting = "cox",
             times = c(36,60),ROC = TRUE)
ROC.test<-timeROC(T=testdata$PFS,delta = testdata$PFS_event,
             marker = testdata$lp,cause = 1,weighting = "cox",
             times = c(36,60),ROC = TRUE)

plot(ROC,time = 36,title = FALSE,lwd=2)
plot(ROC,time = 60,col = 'blue',add = TRUE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 3 years:',round(ROC$AUC[1],2)),
         paste0('AUC at 5 years:',round(ROC$AUC[2],2))),
       col = c('red','blue'),lwd = 2,bty = 'n')

plot(ROC.test,time = 36,title = FALSE,lwd=2)
plot(ROC.test,time = 60,col = 'blue',add = TRUE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 3 years:',round(ROC.test$AUC[1],2)),
         paste0('AUC at 5 years:',round(ROC.test$AUC[2],2))),
       col = c('red','blue'),lwd = 2,bty = 'n')