library(flam)
exp.data<- read.table(file="Expression_counts.txt",header=T)
#now let's calcuate the variance at each locus and compute the coefficient of variation
#pick the top 2000 loci with the largest CV
#first we will only examine loci with all positive expression data
#Variants of this filter can be tried like five zeros allowed if all in one selection regime
pos.vec<- sapply(1:length(exp.data[,1]),function(x) prod(exp.data[x,2:11]>0,na.rm=T))
exp.data.pos<- exp.data[pos.vec,]
cv.data<- sapply(1:length(exp.data.pos[,1]),function(x){
  sd.row<-sd(exp.data.pos[x,2:11])
  mean.row<- mean(as.numeric(exp.data.pos[x,2:11]))
  sd.row/mean.row
})
cv.order<- order(cv.data,decreasing=T)
#cv.order has the location, in exp.data.pos, of the largest, next largest etc cv's.
#So now we will take the largest 2000
exp.data.2000<- exp.data.pos[cv.order[1:2000],]
#To start we will use mortality in cages from day 20-22 from egg.
#this corresponds to adult days 9-11
mort<- c(0.013154555,0.052245068,0.044368497,0.047775774,0.080684659,0.003042255,0.007140445,0.004185332,0.00463155,0.003414137)
#Next make the feature matrix the transpose of exp.data.2000 and use column 1 as the column names
x.exp<- t(exp.data.2000[,2:11])
colnames(x.exp)<- exp.data.2000[,1]
#run flamCV at six different values of alpha, 0.5, 0.6, 0.7, 0.8, 0.9, 1.
# set seed so the same fols are used, start with 5 fold (10 fold gave error)
alpha5.flamcv<- flamCV(x.exp,mort,alpha=0.5,n.fold=5,seed=1)
alpha6.flamcv<- flamCV(x.exp,mort,alpha=0.6,n.fold=5,seed=1)
alpha7.flamcv<- flamCV(x.exp,mort,alpha=0.7,n.fold=5,seed=1)
alpha8.flamcv<- flamCV(x.exp,mort,alpha=0.8,n.fold=5,seed=1)
alpha9.flamcv<- flamCV(x.exp,mort,alpha=0.9,n.fold=5,seed=1)
alpha10.flamcv<- flamCV(x.exp,mort,alpha=1,n.fold=5,seed=1)
#So lets look fits with the log of mort
alpha9log.flamcv<- flamCV(x.exp,log(mort),alpha=0.9,n.fold=5,seed=1)
alpha8log.flamcv<- flamCV(x.exp,log(mort),alpha=0.8,n.fold=5,seed=1)
alpha7log.flamcv<- flamCV(x.exp,log(mort),alpha=0.7,n.fold=5,seed=1)





