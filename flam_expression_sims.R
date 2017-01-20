# Flam simulation study
# This program will simulation expression data for 10 populations
# Indendently expression data for functional loci will be generated and used
# to simulate phenotypes. These data will mimic ACO/CO adult phenotype and
# expression data that exist. From this simulated data FLAM will be run to
# identify the causative loci under a range of FLAM parameter settings
# this is the code used for sim1.
library(flam)
library(foreach)
library(doParallel)
registerDoParallel(cores=10)
ns<- nm<- nw <- 10 #These are the number if strong, moderate and weak loci
Tpheno<- 5 #The number of phenotypes per selection treatment
a1<- 0.5 #The relative effect of moderate loci
a2<- 0.1 #The relative effect of weak loci
Nloci<- 2000 #number of expression loci in the entire data set
#Next we will start a loop to compute Nsim independent data bases
Nsim<- 100
foreach(i=1:Nsim) %dopar% {
  #Next calculate the expression data for the five ACO population at the functional loci
  aco.e.data<- sapply(1:Tpheno,function(x)
    exp(rt((ns+nm+nw),5)*(1/3.5))) #the mean is about 1.05
  aco.p.data<- sapply(1:Tpheno, function(x) {
    ((sum(aco.e.data[1:ns,x])+a1*sum(aco.e.data[(ns+1):(ns+nm),x])+
    a2*sum(aco.e.data[(ns+nm+1):(ns+nm+nw),x]))/(ns+a1*nm+a2*nw)+rnorm(1,0,0.05))})
                                                                #the last rnorm is environmental variance
  #We calculate the CO phenotypes next. The difference is that the expression is
  # shifted up by 1 unit in the CO's.
  co.e.data<- sapply(1:Tpheno,function(x)
    (exp(rt((ns+nm+nw),5)*(1/3.5))+1)) #the mean should be about 2.05
  co.p.data<- sapply(1:Tpheno, function(x) {
    ((sum(co.e.data[1:ns,x])+a1*sum(co.e.data[(ns+1):(ns+nm),x])+
        a2*sum(co.e.data[(ns+nm+1):(ns+nm+nw),x]))/(ns+a1*nm+a2*nw)+rnorm(1,0,0.05))})
  # Next we calculate the genomic expression data.
  # The first ns+nm+nw loci are the functional loci followed by the neutral loci
  # We draw new samples since the expression data is measured on different flies then we used for phenotyping
  aco.exp<- sapply(1:Tpheno,function(x)
    exp(rt(Nloci,5)*(1/3.5)))
  co.exp<- sapply(1:Tpheno,function(x){
   c((exp(rt((ns+nm+nw),5)*(1/3.5))+1), exp(rt((Nloci-ns-nm-nw),5)*(1/3.5)))
  })
  #Now let's create a single matrix with phenotype, loci #, and expression data
  all.data<- cbind(aco.exp,co.exp)
  all.data<- rbind(c(aco.p.data,co.p.data),all.data)
  #Now let's save these results on disk for later analysis
  write(t(all.data),file=paste("sim",1,"rep",i,sep="_"),ncolumns=10)
}

#Next we need to open each of the data files run flam, save results
library(flam)
library(foreach)
library(doParallel)
registerDoParallel(cores=20)
Nsim=100
sim.out<- foreach(i.sim=1:Nsim) %dopar% {
  loci.incl<- matrix(0,nrow=500,ncol=8)
  a.cv<- matrix(0,nrow=2,ncol=10)
  sim.data<- read.table(file=paste("sim",1,"rep",i.sim,sep="_"))
  sim.pheno<- sim.data[1,] #The first row is the phenotype data
  sim.pheno<- as.numeric(sim.pheno)
  sim.gen<- sim.data[-1,] #The genetic data is in all the other rows
  sim.gen<- t(sim.gen)
  colnames(sim.gen)<- 1:2000
  a<- seq(from=0.55,to=1.0,by=0.05)
  for(j.sim in a){
    sim.flamCV<- flamCV(sim.gen,sim.pheno,alpha=j.sim,n.fold=5,seed=1)#we need to use the same seed for each a
    cv.min<- min(sim.flamCV$mean.cv.error)
    cv.error.min<- sim.flamCV$se.cv.error[which(order(sim.flamCV$mean.cv.error,decreasing=TRUE)==50)]
    a.cv[,which(a==j.sim)]<-c(j.sim,cv.min) 
    #count.flam enumerates the frequency that each locus appears in the non-sparse list among the 50 lambda values
    count.flam<- matrix(0,nrow=1,ncol=2000)
    for (i in 1:2000) {
      total.count<- 0
      for (j in 1:50){
        temp<- sim.flamCV$flam.out$non.sparse.list[[j]]
        total.count<- total.count+length(temp[temp==i])
      }
      count.flam[i]<- total.count
    }
        #Now iterate over the 50 lamda values
    sapply(1:50,function(k.sim){
      #Now check the number of strong loci matched
      sim.loci<- sim.flamCV$flam.out$non.sparse.list[[k.sim]]
      if (length(sim.loci)>0){
        strong.loci<- sum(sapply(1:length(sim.loci),function(i) if(any(1:10==sim.loci[i])) 1 else 0))
        mod.loci<- sum(sapply(1:length(sim.loci),function(i) if(any(11:20==sim.loci[i])) 1 else 0))
        weak.loci<- sum(sapply(1:length(sim.loci),function(i) if(any(21:30==sim.loci[i])) 1 else 0))
        neutral.loci<- sum(sapply(1:length(sim.loci),function(i) if(any(31:2000==sim.loci[i])) 1 else 0))} else {
        strong.loci<- 0
        mod.loci<- 0
        weak.loci<- 0
        neutral.loci<- 0}
      temp<- count.flam[1:30]
      temp1<- length(temp[temp>=k.sim])
      temp<- count.flam[31:2000]
      temp2<- length(temp[temp>=k.sim])
      row.num<- (which(a==j.sim)-1)*50+k.sim
      loci.incl[row.num,]<<- c(j.sim,(sim.flamCV$mean.cv.error[k.sim]-cv.min)/cv.error.min,strong.loci,mod.loci,weak.loci,neutral.loci,temp1,temp2)
      }#end of k.sim
      )
  }# end of j.sim
  list(a.cv,loci.incl)
}# end of i.sim

#The next section will analyze and graph the results
#All results are in sim.out[[i]], i=1,,Nsim
#sim.out[[1]][[1]] is the a.cv matrix, 2 x 10
#sim.out[[1]][[2]] is the loci.incl matrix, 500 x 4

#total.vs.a will compute the sum of strong, moderate and weak loci included in the model vs a values
a<- seq(from=0.55,to=1.0,by=0.05)
total.vs.a<- sapply(1:100,function(x){
  temp<- sim.out[[x]][[2]]
  ave<- NULL
  for (i in a){
    temp2<- sapply(1:500,function(j)
      if (temp[j,1]==i)
        temp[j,3]+temp[j,4]+temp[j,5] else NA)
    ave<- c(ave,mean(temp2,na.rm=TRUE))
  }
  ave
})
plot(a,apply(total.vs.a,1,mean),ylim=c(0,20))
# Now calculate the number of neutral loci vs a and store in in neutral.vs.a
neutral.vs.a<- sapply(1:100,function(x){
  temp<- sim.out[[x]][[2]]
  ave<- NULL
  for (i in a){
    temp2<- sapply(1:500,function(j)
      if (temp[j,1]==i)
        temp[j,6] else NA)
    ave<- c(ave,mean(temp2,na.rm=TRUE))
  }
  ave
})
lines(a,apply(neutral.vs.a,1,mean))

#Now for a single value of a let's look at the sum of all loci vs. std difference
total.vs.diff.95<- NULL
for (i in 1:100) {
  temp<- sim.out[[i]][[2]]
  temp<- temp[temp[,1]==a[9],2:5] #a[9]=0.95?
  temp<- cbind(temp[,1],temp[,2]+temp[,3]+temp[,4])
  total.vs.diff.95<- rbind(total.vs.diff.95,temp)
}
#now we can plot the mean standardized diff for each number of loci included
mean.diff.95<-tapply(total.vs.diff.95[,1],total.vs.diff.95[,2],mean)
plot(0:24,mean.diff.95,ylab="Standardized Difference",xlab="Total Included Loci")
#plot a histogram of only those runs that gave 20 or more loci
total.20.95<- total.vs.diff.95[total.vs.diff.95[,2]>19,1]
temp.breaks<- c(0,1,2,3,4,5,6,7,8,9,10.20,30,40,50,140)
hist(total.20.95,breaks=temp.breaks,xlab="Standardized Difference")

#Below  we will isolate each value of "a" and then compute the average distribution of included loci 
#vs number of lambda values.
ave.count<- sapply(1:10,function(x){
  row.start<- (x-1)*50
  temp<- matrix(0,nrow=50,ncol=2)
  dummy<- sapply(1:50,function(y){
    temp2<- 0
    temp3<- 0
    for (i in 1:100){
      temp2<- temp2+sim.out[[i]][[2]][row.start+y,7] #selected loci
      temp3<- temp3+sim.out[[i]][[2]][row.start+y,8] #neutral loci
    }
    temp[y,]<<-c(temp2/100,temp3/100)
  })
  list(temp)
})

ai<- 10
plot(1:50, ave.count[[ai]][,1],xlab="Number of Lambda Values",ylab="Number of Loci Included in x Lambda Values")
lines(1:50,ave.count[[ai]][,2])

#Next for each value of a.temp we calculate the average number of strong, moderate and weak loci included at each lambda value
a.temp<- 0.55 #change for each value of a
temp2<- NULL
for (i in 1:100){
  temp<- sim.out[[i]][[2]]
  temp<- temp[temp[,1]==a.temp,3:5]
  temp2<<- rbind(temp2,temp)
}
#Calculate the mean for each of 50 lambda's
mean.temp<- mean.sel.55<- sapply(1:50,function(i){ #change for each value of a
  temp.s<- 0
  temp.m<- 0
  temp.w<- 0
  for (j in 1:100){
    temp.s<- temp.s + temp2[(j-1)*50+i,1]
    temp.m<- temp.m + temp2[(j-1)*50+i,2]
    temp.w<- temp.w + temp2[(j-1)*50+i,3]
  }
  c(temp.s/100,temp.m/100,temp.w/100)
})
plot(1:50,mean.temp[1,],type="l",lty=1,ylab="Number of included loci (out of 10)",xlab="lamda sequence")
lines(1:50,mean.temp[2,],lty=2)
lines(1:50,mean.temp[3,],lty=3)
