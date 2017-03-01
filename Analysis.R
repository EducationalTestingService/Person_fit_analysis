irt=function(t,a,b){return(1/(1+exp(a*(b-t))))}
# Simulating the true ability for an aberrant examinee
aberth=function(n,ab){switch(ab,ch=-2+1.5*runif(n),cr=0.5+1.5*runif(n),g=-2+1.5*runif(n),ca=0.5+1.5*runif(n),r=-2+4*runif(n))}
#Function to compute response probabilities for an aberrant examinee
aber=function(a,b,t,ab){switch(ab,ch=ifelse(b<1.5,irt(t,a,b),1),cr=ifelse(b<= -1.5,0,irt(t,a,b)),g=ifelse(b<0.5,irt(t,a,b),0.25),ca=ifelse(b<= -0.5,0.5,irt(t,a,b)),r=rep(0.25,length(b)))}
#R packages MESS, metafor, PerFit,ltm and irtoys should be installed
library(MESS)
library(metafor)
library(PerFit)
library(ltm)
library(irtoys)
#Function to compute an ROC Area and standard error
ROC=function(good,bad,x) 
{ Areas=rep(0,ncol(good))
  SE=Areas  
  F=matrix(0,length(x),ncol(good))
  H=F
  for (j in 1:ncol(good))
  {goodstat=good[,j]
   badstat=bad[,j]
    for (i in 1:length(x)){
  F[i,j]=length(goodstat[goodstat<x[i]])/length(goodstat)
   H[i,j]=length(badstat[badstat<x[i]])/length(badstat)}
   #Use the R function `auc' to compute the ROC area
   a=auc(F[,j],H[,j],type='spline')
   Areas[j]=ifelse(a>1,0.999,a)
   a=Areas[j]
   ng=nrow(good)
   nb=nrow(bad)
   q1=a/(2-a)
   q2=2*a*a/(1+a)
   SE[j]=sqrt(( a*(1-a)+(nb-1)*(q1-a*a)+(ng-1)*(q2-a*a) )/(ng*nb))}
   return(rbind(Areas,SE,F,H))}
Output=NULL
SEs=NULL
x=seq(-5,5,,100)
nexam=10000
types=c("ch","cr","g","ca","r")
statsK03=NULL#PFSs of all examinees over 60 data sets
aberall=NULL#Indicator of aberrance of all examinees
areaitems=NULL#ROC areas for different test lengths
cat("#item , type,   %Aber,   lz,   HT,   U3 \n",file="RAreas_PnP",append=FALSE)
for (nitem in c(17,33,65)){
yitems=NULL#Item scores of all examinees for a given test length
aberitems=NULL#Indicator of aberrance of all examinees for a given test length
 for (nty in 1:5){
for (propaber in c(0.05,0.1,0.25,0.5)){
Exam=1:nexam
type=types[nty]    
y=matrix(0,nexam,nitem)#Scores of all examinees for a test length,%Aber,Aber-type 
naber=nexam*propaber#Number of aberrant examinees in a data set
theta=runif(nexam,-2,2)#Simulate true ability of non-aberrant examinees
theta[1:naber]=aberth(naber,type)#True ability of aberrant examinees
at=rep(1,nitem)
bt=-2+4*(1:nitem - 1)/(nitem-1)
# Simulate the scores of the examinees
for (i in 1:naber){p = aber(at,bt,theta[i],type)#Aberrant examinees
y[i,]=rbinom(nitem,size=1,prob=p)}
for (i in (naber+1):nexam){p=irt(theta[i],at,bt)#Non-aberrant examinees
     y[i,]=rbinom(nitem,size=1,prob=p)}
scores.rasch = rasch(y) #Fit the Rasch model using package 'ltm'
itparm=coef(scores.rasch)#Estimated item parameters
b=itparm[,1]#Difficulty parameters
a=itparm[,2]#Common slope parameter
raw=apply(y,1,sum)
y=y[raw>0 & raw<nitem,]#Remove those with 0 or full score: HT undefined for them
theta=theta[raw>0 & raw<nitem]
Exam=Exam[raw>0 & raw<nitem]
yitems=rbind(yitems,y)
abe = rep(1,nrow(y))
abe[Exam>naber]=0
aberitems=c(aberitems,abe)
aberall=c(aberall,abe)
#Estimate the examinee abilities using package 'irtoys'
itparms=cbind(a,b,rep(0,nitem))
thetaest=mlebme(y,itparms)[,1]
#Compute lz, HT, and U3 using R package Perfit
lzs=lz(y,Ability=thetaest,IP=itparms)
Hts = Ht(y)
U3s=U3(y,Ability=thetaest,IP=itparms)#U3 is large for aberrant examinees
stats=cbind(lzs[[1]]$PFscores,Hts[[1]]$PFscores,-U3s[[1]]$PFscores)
statsK03=rbind(statsK03,stats)
stats=scale(stats,center=TRUE,scale=TRUE)#Standardize each PFS
good=stats[Exam>naber,]#PFSs for the non-aberrant examinees
bad=stats[Exam<=naber,]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:(length(x)+2),]#False-alarm rates
H=ROCOut[(length(x)+3):nrow(ROCOut),]#Hit rates
SEs=rbind(SEs,ROCOut[2,])
Output=rbind(Output,c(nitem,nty,100*propaber,Areas))
}}
# Suggested approach to compute ROC Area for each test length aggregating
# over the percentages and types of aberrance
scores.rasch = rasch(yitems) #Fit the Rasch model to the combined data
itparm=coef(scores.rasch)#Estimated item parameters
b=itparm[,1]#Difficulty parameters
a=itparm[,2]#Common slope parameter
itparms=cbind(a,b,rep(0,nitem)) 
#Compute the ability estimates for the combined data using `irtoys'
thetaest=mlebme(yitems,itparms)[,1]
Hts = Ht(yitems)
lzs=lz(yitems,Ability=thetaest,IP=itparms)
U3s=U3(yitems,Ability=thetaest,IP=itparms)
stats=cbind(lzs[[1]]$PFscores,Hts[[1]]$PFscores,-U3s[[1]]$PFscores)
stats=scale(stats,center=TRUE,scale=TRUE)
good=stats[aberitems==0,]#PFSs of all nonaberrant examinees for a test length
bad=stats[aberitems==1,]#PFSs of all aberrant examinees for a test length
ROCOut=ROC(good,bad,x)
Areas=ROCOut[1,]#ROC areas for a test length
areaitems=rbind(areaitems,Areas)}
write.table(round(Output,2),"AUC",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)#Write the ROC area for the 60 data sets in file `AUC'
Areas=Output[,4:6]
AUCs=NULL
for (i in 1:ncol(Areas))#Use DerSimonian-Laird algorithm
     {AUCs=c(AUCs,rma(Areas[,i],SEs[,i],method="DL")$b)}
 cat("\n Overall ROC areas Using DerSimonian-Laird algorithm :",round(AUCs,2),"\n",file="AUC",append=TRUE)
good=statsK03[aberall==0,]#PFSs of all nonaberrant examinees 
bad=statsK03[aberall==1,]#PFSs of all aberrant examinees 
ROCOut=ROC(good,bad,x)
cat("\n Overall ROC areas Using Karabatsos (2003) Approach :",round(ROCOut[1,],2),"\n",file="AUC",append=TRUE)
cat("\n\n Overall ROC Areas using the suggested approach for each length:","\n",file="AUC",append=TRUE)
cat("#items , lz,     HT,     U3 \n",file="AUC",append=TRUE)
write.table(cbind(c(17,33,65),round(areaitems,2)),"AUC",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
