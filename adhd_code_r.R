

# The code is suitable for our paper "Abnormal hemispheric asymmetry of
# both brain function and structure in attention deficit/hyperactivity
# disorder: a meta-analysis of individual participant data".
# 2021-01-15

setwd("F:/asymmestry_all/")
gmv7p=read.csv('F:/asymmestry_all/gmv7p.csv',header = FALSE)
gmv7ptoz=matrix(data = 0, nrow = 56, ncol = 7)
for (i in 1:56) {
  for (k in 1:7) {
    gmv7ptoz[i,k]=qnorm(gmv7p[i,k]/2, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    
  }
}
write.csv(gmv7ptoz,'gmv7ptoz.csv')
#
gmv7ptoz2=read.csv('F:/asymmestry_all/gmv7ptoz2.csv',header = FALSE)
n=matrix(data=0,nrow=7,ncol=1)
n[1,1]=47
n[2,1]=66
n[3,1]=241
n[4,1]=78
n[5,1]=90
n[6,1]=64
n[7,1]=41
fenzi=matrix(data=0,nrow=56,ncol=1)
fenmu=matrix(data=sqrt(627),nrow=56,ncol=1)
for (i in 1:56) {
  for (k in 1:7) {
    fenzi[i,1]=fenzi[i,1]+sqrt(n[k,1])*gmv7ptoz2[i,k]
  }
}  
combineZ=fenzi/fenmu
combineZtoP=matrix(data=0,nrow=56,ncol=1)
for (i in 1:56) {
  combineZtoP[i,1]=2*pnorm(abs(combineZ[i,1]), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)#   pval=2P{A>|z|}
}
write.csv(cbind(combineZ,combineZtoP),'combineZPval.csv')
padj=p.adjust(combineZtoP, method = "fdr", n = length(combineZtoP))
sn=which(padj<0.05)
snsn=matrix(sn)
write.csv(snsn,'sig_roi.csv')
#
Q=read.csv('F:/asymmestry_all/Q.csv',header = FALSE)
QtoP=matrix(data=0,nrow=56,ncol=1)
for (i in 1:123) {
  QtoP[i,1]=pchisq(Q[i,1], df=6,ncp=0, lower.tail = FALSE, log.p = FALSE)
}
write.csv(cbind(Q,QtoP),'QPval.csv')

