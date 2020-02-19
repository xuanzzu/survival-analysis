setwd("I://脊索瘤预后预测模型//结果//家琦分析结果//随机筛选指标结果//手术入路 + 斜坡部位 + degree_of_resection    蛋白：e_cad_cyto + Ki_67_2_n + VEGFA_cyto VEGFAp值更小//PFS")
#读取数据
traindata<-read.table("trainData.csv",header = T,sep=",")
library(survival)
library(rms)
library(plyr)
#1 建立生存对象
y<-Surv(traindata$PFS,traindata$PFS_event==1)
#2 单因素COX回归
#2.1 每次只做一个因素的COX回归
unicox<-coxph(y~年龄,data = traindata)
unicox.sum<-summary(unicox)
CI<-paste0(round(unicox.sum$conf.int[,3:4],2),collapse = "-")
HR<-round(unicox.sum$coefficients[,2],2)
Pvalue<-round(unicox.sum$coefficients[,5],4)
unicox.result<-data.frame("characteristics"="Age",
                          "Hazard Ratio"=HR,
                          "CI95"=CI,
                          "P value"=Pvalue)
#2.2 批量做每个因素的COX回归
Unicox.result<-function(x){
  FML<-as.formula(paste0("y~",x))
  unicox<-coxph(FML,ties=c("breslow"),data = traindata)
  unicox.sum<-summary(unicox)
  CI<-paste0(round(unicox.sum$conf.int[,3:4],2),collapse = "-")
  HR<-round(unicox.sum$coefficients[,2],2)
  Pvalue<-round(unicox.sum$coefficients[,5],4)
  unicox.result<-data.frame("characteristics"=x,
                            "Hazard Ratio"=HR,
                            "CI95"=CI,
                            "P value"=Pvalue)
  return(unicox.result)
}
VarNames<-colnames(traindata)[c(3:18,24:83)]
Univar<-lapply(VarNames,Unicox.result)
Univar<-ldply(Univar,data.frame)
#3 多因素COX回归
#3.1 筛选单因素p小于0.05的变量
f<-Univar$characteristics[Univar$P.value<0.05]
fml1<-as.formula(paste0("y~",paste0(f[c(3,4,6)],collapse = "+")))
fml2<-as.formula(paste0("y~",paste0(f[8:15],collapse = "+")))
#3.2 多因素COX回归
MultiCox<-coxph(fml1,ties=c("breslow"),data = traindata)
MultiStep<-step(MultiCox,direction = "backward")
MultiSum<-summary(step(MultiCox,direction = "backward"))

MultiCox.biomarker<-coxph(fml2,ties=c("breslow"),data = traindata)
MultiStep.biomarker<-step(MultiCox.biomarker,direction = "backward")

#3.3 最终模型
MultiCox.final<-coxph(y~ degree_of_resection + e_cad_cyto + Ki_67_2_n + VEGFA_cyto,ties = c("breslow"),
                      data = traindata)

#4 绘制列线图
#4.1 多因素COX回归（绘制列线图用）
MultiCox.nomogram<-cph(y~degree_of_resection + e_cad_cyto + Ki_67_2_n + VEGFA_cyto,method = c("breslow"),
                       data = traindata,x=T,y=T,surv=T)
MultiCox2<-cph(y~degree_of_resection + e_cad_cyto + Ki_67_2_n + VEGFA_cyto,data=traindata,x=T,y=T,surv=T)
MultiStep2<-step(MultiCox2,direction = "backward")
MultiSum2<-summary(step(MultiCox2,direction = "backward"))
#4.2 数据打包
dd<-datadist(traindata)
options(datadist = "dd")
#4.3 生成函数
surv<-Survival(MultiCox.nomogram)
#4.4 建立nomogram
surv1<-function(x)surv(3*12,lp=x)
surv2<-function(x)surv(5*12,lp=x)
nom<-nomogram(MultiCox.nomogram,fun=list(surv1,surv2),lp=F,
              funlabel = c("3-year progression-free survival","5-year progression-free survival"),
              maxscale = 100,
              fun.at = c(0.95,0.9,0.8,0.6,0.4,0.2,0.05))
plot(nom)

#5 内部验证
#5.1 区分度（discrimination）。计算C指数（C-index=1-C）
validate(MultiCox.nomogram,method = "boot",B=1000,dxy = T)
rcorrcens(y~predict(MultiCox.nomogram),data = traindata)
#5.2 一致性（calibration）
#5.2.1 3-year progression-free survival calibartion curve
f3<-cph(y~degree_of_resection + e_cad_cyto + Ki_67_2_n + VEGFA_cyto,method = c("breslow"),
        data = traindata,x=T,y=T,surv=T,time.inc = 36)
cal3<-calibrate(f3,cmethod = "KM",method = "boot",u=36,m=50,B=1000)
plot(cal3,xlab = "Predicted 36 months survival", ylab = "Fraction surviving 36 months")
#5.2.2 5-year progression-free survival calibartion curve
f5<-cph(y~degree_of_resection + e_cad_cyto + Ki_67_2_n + VEGFA_cyto,method = c("breslow"),
        data = traindata,x=T,y=T,surv=T,time.inc = 60)
cal5<-calibrate(f5,cmethod = "KM",method = "boot",u=60,m=50,B=1000)
plot(cal5,xlab = "Predicted 60 months survival", ylab = "Fraction surviving 60 months")

#6 外部验证
testdata<-read.table("testData.csv",header = T,sep=",")
#6.1 区分度（discrimination）
ftest<-cph(Surv(PFS,PFS_event)~predict(MultiCox.nomogram,newdata = testdata),method = c("breslow"),
           x=T,y=T,surv=T,data = testdata)
validate(ftest,method = "boot",B=1000,dxy = T)
rcorrcens(Surv(PFS,PFS_event)~predict(MultiCox.nomogram,newdata = testdata),
          data = testdata)
#6.2 一致性（calibration）
#6.2.1 3-year overall survival calibartion curve
ftest3<-cph(Surv(PFS,PFS_event)~predict(MultiCox.nomogram,newdata = testdata),method = c("breslow"),
            x=T,y=T,surv=T,data = testdata,time.inc = 36)
caltest3<-calibrate(ftest3,cmethod = "KM",method = "boot",u=36,m=16,B=1000)
plot(caltest3,xlab = "Predicted 36 months survival", ylab = "Fraction surviving 36 months")
#6.2.2 5-year overall survival calibartion curve
ftest5<-cph(Surv(PFS,PFS_event)~predict(MultiCox.nomogram,newdata = testdata),method = c("breslow"),
            x=T,y=T,surv=T,data = testdata,time.inc = 60)
caltest5<-calibrate(ftest5,cmethod = "KM",method = "boot",u=60,m=16,B=1000)
plot(caltest5,xlab = "Predicted 60 months survival", ylab = "Fraction surviving 60 months")

#7 计算每个患者评分
Total_Points<-(traindata$degree_of_resection-1)*28.1512+(traindata$e_cad_cyto-60)*0.9090909+traindata$Ki_67_2_n*1.396713+(traindata$VEGFA_cyto-40)*0.1091895
traindata$Total_Points<-Total_Points

Total_Points_testdata<-(testdata$degree_of_resection-1)*28.1512+(testdata$e_cad_cyto-60)*0.9090909+testdata$Ki_67_2_n*1.396713+(testdata$VEGFA_cyto-40)*0.1091895
testdata$Total_Points<-Total_Points_testdata
#7.1 对比低风险组和高风险组PFS
traindata_with_points<-read.table("traindata with points.csv",header = T,sep=",")
KM<-survfit(Surv(PFS,PFS_event)~Total_Points二分类,data=traindata_with_points)
KMdiff<-survdiff(Surv(PFS,PFS_event==1)~Total_Points二分类,data=traindata_with_points)

testdata_with_points<-read.table("testdata with points.csv",header = T,sep=",")
KMtest<-survfit(Surv(PFS,PFS_event)~Total_Points二分类,data=testdata_with_points)
KMdiff<-survdiff(Surv(PFS,PFS_event==1)~Total_Points二分类,data=testdata_with_points)
#7.2 绘制K-M曲线（low-risk group vs high-risk group)
library(survminer)
ggsurvplot(KM,legend.labs=c("low-risk group","high-risk group"), data = traindata_with_points)
ggsurvplot(KMtest,legend.labs=c("low-risk group","high-risk group"), data = testdata_with_points)

#Table 3的计算
t.test(Total_Points~Total_Points二分类,data=traindata_with_points)
sd(traindata_with_points$Total_Points[which(traindata_with_points$Total_Points二分类==1)])
