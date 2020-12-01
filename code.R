#code for paper: Accuracy Analysis of Choropleth Maps: A Perspective from Attribute Uncertainty
#（1）unifrom distributed samples
MCJY<-function(x,m,ra) 
{
l<-length(x)
y<-1:(m*l)
dim(y)<-c(l,m)
for(i in 1:l)y[i,]<-runif(m,x[i]-ra*x[i],x[i]+ra*x[i])
for(j in 1:m)y[,j]<-sort(y[,j])
write.csv(y,file="exprimental_data_for_section3\\uniform_distribution_124.csv",row.names=FALSE)
}

#（2）normal distributed samples
MCZHT<-function(x,m,ra) 
{
l<-length(x)
sigma<-diag(1,l,l)
for(i in 1:l)
{
sigma[i,i]<-(x[i]*ra)* (x[i]*ra)
}
y<-mvrnorm(n=m,x,sigma)
for(i in 1:m)y[i,]<-sort(y[i,])
z<-t(y)
write.csv(z,file="exprimental_data_for_section3\\normal_distribution_124.csv",row.names=FALSE)
}

#（3）Equal Intervals
EIGROUP<-function(x,n) 
{
group <- vector("list", n)
for(i in 1:n)group[[i]]<-x[cut(x,n,labels=F)==i]
group
}

#（4）Quantile
QGROUP<-function(x,n)
{
m<-length(x)
s<-floor(m/n)
t<-m%%n
group <- vector("list", n)
if(t==0)
{
for(i in 1:n)group[[i]]<-x[((i-1)*s+1):(i*s)]
}
else
{
group[[1]]<-x[1:s]
o<-x[(s+1):((t+1)*s+t)]
for(i in 1:t)group[[i+1]]<-o[((i-1)*(s+1)+1):(i*(s+1))]
}
if(t+1<n)
{
p<-x[((t+1)*s+t+1):m]
for(i in 1:(n-t-1))group[[i+t+1]]<-p[((i-1)*s+1):(i*s)]
}
group
}

#（5）Jenks）（Natural Breaks）
NBGROUP<-function(x,n)
{
group<-QGROUP(x,n)
group1<-NB1group(x,group,n)
group2<-forcegroup(x,group1,n)
group2
}
diedaigroup<-function(x,group,n)
{l<-length(x)
m.group<-c()
nb.group<- vector("list", n)
d.mean<-c(rep(0,l*n))
dim(d.mean)<-c(l,n)
for(i in 1:n)
{
m.group[i]<-mean(group[[i]])
d.mean[,i]<-abs(x-m.group[i])
}
for(i in 1:l)
{
xb<-which(d.mean[i,]==min(d.mean[i,]))
nb.group[[xb]]<-c(nb.group[[xb]],x[i])
}
nb.group
}
NB1group<-function(x,group,n)
{
NB1.group<-diedaigroup(x,group,n)
tai1<-taigroup(x,group,n)
tai2<-taigroup(x,NB1.group,n)
while(tai2>tai1)
{group<-NB1.group
NB1.group<-diedaigroup(x,group,n)
tai1<-taigroup(x,group,n)
tai2<-taigroup(x,NB1.group,n)}
group
}
forcegroup1<-function(x,group,n)
{
f.group<-group
for(i in 1:(n-1))
{
tai1<-taigroup(x,group,n)
f.group[[i]]<-c(f.group[[i]],f.group[[i+1]][1])
u<-length(f.group[[i+1]])
f.group[[i+1]]<-f.group[[i+1]][2:u]
tai2<-taigroup(x,f.group,n)
while(tai2>tai1)
{group<-f.group
tai1<-taigroup(x,group,n)
f.group[[i]]<-c(f.group[[i]],f.group[[i+1]][1])
u<-length(f.group[[i+1]])
f.group[[i+1]]<-f.group[[i+1]][2:u]
tai2<-taigroup(x,f.group,n)}
f.group<-group
}
group
}
forcegroup2<-function(x,group,n)
{
f.group<-group
for(i in 0:(n-2))
{
tai1<-taigroup(x,group,n)
u<-length(f.group[[n-i-1]])
f.group[[n-i]]<-c(f.group[[n-i-1]][u],f.group[[n-i]])
f.group[[n-i-1]]<-f.group[[n-i-1]][1:(u-1)]
tai2<-taigroup(x,f.group,n)
while(tai2>tai1)
{group<-f.group
tai1<-taigroup(x,group,n)
u<-length(f.group[[n-i-1]])
f.group[[n-i]]<-c(f.group[[n-i-1]][u],f.group[[n-i]])
f.group[[n-i-1]]<-f.group[[n-i-1]][1:(u-1)]
tai2<-taigroup(x,f.group,n)}
f.group<-group
}
group
}
forcegroup<-function(x,group,n)
{
group1<-forcegroup1(x,group,n)
group2<-forcegroup2(x,group1,n)
group3<-forcegroup2(x,group2,n)
tai2<-taigroup(x,group2,n)
tai3<-taigroup(x,group3,n)
while(tai3>tai2)
{group2<-group3
group3<-forcegroup2(x,group2,n)
tai2<-taigroup(x,group2,n)
tai3<-taigroup(x,group3,n)}
nb2group<-NB1group(x,group2,n)
nb2group
}

#（6）（TAI）、（GVF）
taigroup<-function(x,group,n)
{
l<-length(x)
m.group<-c()
e.group<- vector("list", n)
cs.e.group<-c()
all.e.group<-c()
for(i in 1:n)
{
m.group[i]<-mean(group[[i]])
s<-length(group[[i]])
j<-1
while(j<s+1)
{e.group[[i]][j]<-abs(group[[i]][j]-m.group[i])
j<-j+1}
cs.e.group[i]<-sum(e.group[[i]])
}
all.e.group<-sum(cs.e.group)
z<-abs(x-mean(x))
tai<-1- all.e.group/sum(z)
tai
}

gvfgroup<-function(x,group,n)
{
l<-length(x)
m.group<-c()
e.group<- vector("list", n)
cs.e.group<-c()
all.e.group<-c()
for(i in 1:n)
{
m.group[i]<-mean(group[[i]])
s<-length(group[[i]])
j<-1
while(j<s+1)
{e.group[[i]][j]<- (group[[i]][j]-m.group[i])*(group[[i]][j]-m.group[i])
j<-j+1}
cs.e.group[i]<-sum(e.group[[i]]) 
}
all.e.group<-sum(cs.e.group)
z<-(x-mean(x))*(x-mean(x))
gvf<-1- all.e.group/sum(z)
gvf
}

#（7）（GVF）
yangbenEIgvf<-function(x,tg) 
{
m<-length(x[1,])
yangbentai<-matrix(rep(0, m* (tg-1)),nrow=tg-1,ncol=m)
danyangbentai<-matrix(rep(0, (tg-1)))
for(i in 1:m)
{g<-2
while(g<tg+1)
{group<-EIGROUP(x[,i],g)
danyangbentai[g-1]<-gvfgroup(x[,i],group,g)
g<-g+1}
yangbentai[,i]<-danyangbentai}
write.csv(yangbentai,file="C:\\xufeng\\fenjimapshiyan\\yangbenEIgvf.csv",row.names=FALSE)
}

yangbenQgvf<-function(x,tg) 
{
m<-length(x[1,])
yangbentai<-matrix(rep(0, m* (tg-1)),nrow=tg-1,ncol=m)
danyangbentai<-matrix(rep(0, (tg-1)))
for(i in 1:m)
{g<-2
while(g<tg+1)
{group<-QGROUP(x[,i],g)
danyangbentai[g-1]<-gvfgroup(x[,i],group,g)
g<-g+1}
yangbentai[,i]<-danyangbentai}
write.csv(yangbentai,file="C:\\xufeng\\fenjimapshiyan\\yangbenQgvf.csv",row.names=FALSE)
}
yangbenNBgvf<-function(x,tg) 
{
m<-length(x[1,])
yangbentai<-matrix(rep(0, m* (tg-1)),nrow=tg-1,ncol=m)
danyangbentai<-matrix(rep(0, (tg-1)))
for(i in 1:m)
{g<-2
while(g<tg+1)
{group<-NBGROUP(x[,i],g)
danyangbentai[g-1]<-gvfgroup(x[,i],group,g)
g<-g+1}
yangbentai[,i]<-danyangbentai}
write.csv(yangbentai,file="C:\\xufeng\\fenjimapshiyan\\yangbenNBgvf.csv",row.names=FALSE)
}

#output
#（1）Monte Carlo for uniform distribution
x<-as.matrix(read.csv("C:\\xufeng\\fenjimapshiyan\\junyunfenbu124.csv"))
x1<-MCJY(x,1000,2/100) 
y1<-MCZHT(x,1000,2/100) 
x2<-MCJY(x,1000,4/100) 
y2<-MCZHT(x,1000,4/100)
x3<-MCJY(x,1000,6/100) 
y3<-MCZHT(x,1000,6/100)
x4<-MCJY(x,1000,8/100) 
y4<-MCZHT(x,1000,8/100)
x5<-MCJY(x,1000,10/100) 
y5<-MCZHT(x,1000,10/100)
x6<-MCJY(x,1000,12/100) 
y6<-MCZHT(x,1000,12/100)
x7<-MCJY(x,1000,14/100) 
y7<-MCZHT(x,1000,14/100)
x8<-MCJY(x,1000,16/100) 
y8<-MCZHT(x,1000,16/100)
x9<-MCJY(x,1000,18/100) 
y9<-MCZHT(x,1000,18/100)
x10<-MCJY(x,1000,20/100) 
y10<-MCZHT(x,1000,20/100)
#（2）Monte Carlo for normal distribution
x<-as.matrix(read.csv("C:\\xufeng\\fenjimapshiyan\\正态分布数据\\zhengtaifenbu124.csv"))
x1<-MCJY(x,1000,2/100) 
y1<-MCZHT(x,1000,2/100) 
x2<-MCJY(x,1000,4/100) 
y2<-MCZHT(x,1000,4/100)
x3<-MCJY(x,1000,6/100) 
y3<-MCZHT(x,1000,6/100)
x4<-MCJY(x,1000,8/100) 
y4<-MCZHT(x,1000,8/100)
x5<-MCJY(x,1000,10/100) 
y5<-MCZHT(x,1000,10/100)
x6<-MCJY(x,1000,12/100) 
y6<-MCZHT(x,1000,12/100)
x7<-MCJY(x,1000,14/100) 
y7<-MCZHT(x,1000,14/100)
x8<-MCJY(x,1000,16/100) 
y8<-MCZHT(x,1000,16/100)
x9<-MCJY(x,1000,18/100) 
y9<-MCZHT(x,1000,18/100)
x10<-MCJY(x,1000,20/100) 
y10<-MCZHT(x,1000,20/100)
#（3）Monte Carlo for skewed distribution
x<-as.matrix(read.csv("C:\\xufeng\\fenjimapshiyan\\piantaifenbu124.csv"))
x1<-MCJY(x,1000,2/100) 
y1<-MCZHT(x,1000,2/100) 
x2<-MCJY(x,1000,4/100) 
y2<-MCZHT(x,1000,4/100)
x3<-MCJY(x,1000,6/100) 
y3<-MCZHT(x,1000,6/100)
x4<-MCJY(x,1000,8/100) 
y4<-MCZHT(x,1000,8/100)
x5<-MCJY(x,1000,10/100) 
y5<-MCZHT(x,1000,10/100)
x6<-MCJY(x,1000,12/100) 
y6<-MCZHT(x,1000,12/100)
x7<-MCJY(x,1000,14/100) 
y7<-MCZHT(x,1000,14/100)
x8<-MCJY(x,1000,16/100) 
y8<-MCZHT(x,1000,16/100)
x9<-MCJY(x,1000,18/100) 
y9<-MCZHT(x,1000,18/100)
x10<-MCJY(x,1000,20/100) 
y10<-MCZHT(x,1000,20/100)

#（4）output of GFV
x<-as.matrix(read.csv("C:\\xufeng\\fenjimapshiyan\\正态分布数据\\ZHTFBjunyunfenbuyangben20-100.csv"))
EIgvf<-yangbenEIgvf(x,7)
Qgvf<-yangbenQgvf(x,7)
NBgvf<-yangbenNBgvf(x,7)

