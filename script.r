setwd("~/particarlos/pythia8303/rpythia")
source('ler.r')


s=mesmo_proceso('s',0)
h=mesmo_proceso('h',1)
w=mesmo_proceso('w',2)
t=mesmo_proceso('t',3)

save.image(file='myEnvironment_eta_oitocentos.RData')



#library(transport)
library(RcppArmadillo)
library(Rcpp)
# Use RcppXPtrUtils for simple usage of C++ external pointers
library(RcppXPtrUtils)
library(parallelDist)

## user-defined distance function
minkowskiFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
  return (arma::accu(arma::square(A - B))-2*pow(A(0,0)-B(0,0),2));
}", depends = c("RcppArmadillo"))



trivial_dist  <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
  return ((A(0,0)!=B(0,0))+pow(A(0,1)-B(0,1),2)+(A(0,2)!=B(0,2))+(A(0,3)!=B(0,3))+(A(0,4)!=B(0,4)));
}", depends = c("RcppArmadillo"))





# cppFunction('double ip(arma::vec p, arma::vec v) {
#   double impact = arma::norm(arma::cross(p,v));
#   return impact;
# }',depends = c("RcppArmadillo"))

cppFunction('arma::vec ips(arma::mat P, arma::mat V) {
  int n=P.n_rows;
  arma::vec p(3);
  arma::vec v(3);
  arma::vec impacts(n);
  double np=0;
  
  for (int i = 0; i < n; i++) {
  p=P.row(i).t();
  v=V.row(i).t();
  np=arma::norm(p);
  if(np<1e-16) {
  impacts(i)=0;
  }else {
  impacts(i)=arma::norm(arma::cross(p,v)); 
  impacts(i)=impacts(i)/np;
  }
}
  return impacts;
}',depends = c("RcppArmadillo"))

# o parametro de impacto coincide coa magnitude do momento angular no production vertex da neta 
# partido pola norma de p https://inis.iaea.org/collection/NCLCollectionStore/_Public/49/103/49103732.pdf


estables=rbind(s,h,w,t)

names(estables)[names(estables) == "rapi"] <- "eta"

#estables$el=(estables$e)

lista=aggregate(estables$proc, list(estables$ev),max)




#ipss=ips(cbind(estables$p_x,estables$p_y,rep(0,length(estables$p_x))),cbind(estables$x,estables$y,rep(0,length(estables$p_x))))
ipss=ips(cbind(estables$p_x,estables$p_y,estables$p_z),cbind(estables$x,estables$y,estables$z))


estables$ips=ipss
estables$mag_pv=sqrt(estables$x^2+estables$y^2+estables$z^2)
estables_l=estables
estables=estables_l[estables$mag_pv<5e2,]


# estables$n1=estables$e/estables$et
# estables$n2=estables$p_x/estables$et
# estables$n3=estables$p_y/estables$et
# estables$n4=estables$p_z/estables$et

#grupo=as.matrix(cbind(estables$n1,estables$n2,estables$n3,estables$n4))

#sigma=cov(grupo)
#sigmainv=solve(sigma)

#sigmam=chol(sigmainv)
#grupo=grupo%*%t(sigmam)


# estables$n1=grupo[,1]
# estables$n2=grupo[,2]
# estables$n3=grupo[,3]
# estables$n4=grupo[,4]
# normalize=function(x) {
#   return((x-min(x))/(max(x)-min(x)))
# }


normalizacions=aggregate(estables$e, list(estables$ev), sum)
nev=nrow(normalizacions)

#estables$en=estables$el


for(i in 1:nev) {
  estables$en[estables$ev==i]=estables$el[estables$ev==i]/normalizacions$x[i]
}



d=matrix(0,nrow=nev,ncol=nev)
for(k in 1:nev) {
  if(k==nev) {
    break
  }
  for(p in (k+1):nev) {
    
   
    ind1=(estables$ev==p)
    ind2=(estables$ev==k)
    
    cuadri1=estables[ind1,c('pt','rapi','phi')]
    cuadri2=estables[ind2,c('pt','rapi','phi')]
    #cuadri1=apply(cuadri1,2,normalize)
    #cuadri2=apply(cuadri2,2,normalize)

#disc1=estables[ind1,c('q','spin','b','le','lm')]
#disc2=estables[ind2,c('q','spin','b','le','lm')]



ip1=estables[ind1,c('ips')]
ip2=estables[ind2,c('ips')]
#ip1=apply(ip1,2,normalize)
#ip2=apply(ip2,2,normalize)



en1=estables[ind1,c('en')]
en2=estables[ind2,c('en')]



n=nrow(cuadri1)
m=nrow(cuadri2)
custos=matrix(0,nrow=n,ncol=m)


# minkowski<-function(cuadri1,cuadri2) {
#   return(-cuadri1[1]*cuadri1[1]+sum(cuadri1[2:4]*cuadri1[2:4])-cuadri2[1]*cuadri2[1]+sum(cuadri2[2:4]*cuadri2[2:4])-2*(-cuadri1[1]*cuadri2[1]+sum(cuadri1[2:4]*cuadri2[2:4])))
# }
# 
# for(i in 1:n) {
#   for(j in 1:m){
#     
#     custos[i,j]=minkowski(as.numeric(proba1[i,]),as.numeric(proba2[j,]))
# }
# }

res1=parDist(as.matrix(rbind(cuadri1,cuadri2)), method="custom", func = minkowskiFuncPtr)
#res1=parDist(as.matrix(rbind(cuadri1,cuadri2)), method="euclidean")
#res2=parDist(as.matrix(rbind(disc1,disc2)), method="custom", func = trivial_dist)
if((sum(ip1)<1e-20) & (sum(ip2)<1e-20) ){res3=matrix(rep(0,(n+m)^2),ncol=n+m)
}else{res3=parDist(as.matrix(c(ip1,ip2)),method='euclidean')}

castres=sqrt(sqrt(abs(as.matrix(res1)))+(as.matrix(res3)))#+sqrt(as.matrix(res2))



# }





custos=t(castres[-seq(1,n),-seq(n+1,m+n)])


d[k,p]=wasserstein(as.numeric(en1), as.numeric(en2), p=2, tplan=NULL, costm=abs(custos),prob=TRUE)


  }  
}
library(MASS)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

mpt=max(log(estables$pt+1))

estables$lptn=log(estables$pt+1)/mpt
estables$lipsn=log(estables$ips+1)/max(log(estables$ips+1))
estables$etan=estables$eta/max(abs(estables$eta))/sqrt(2)
estables$phin=estables$phi/max(abs(estables$phi))/sqrt(2)

for(k in 201:800){
 
    ind1=(estables$ev==k & estables$lips>0)
    #cuadri1=estables[ind1 & estables$ips>0 & estables$pt>0, c('pt','ips')]
    cuadri1=estables[ind1,c('pt','etan','phin','lipsn')]
    x=sqrt(cuadri1$etan^2+cuadri1$phin^2)
    
    y=cuadri1$lipsn
    
   
    # Default call 
    #ker <- kde2d(log(cuadri1$pt), log(cuadri1$ips),n=20,lims=c(c(-4,4),c(-35,5)))
    ker <- kde2d(x, y,n=15,lims=c(c(0,1),c(0,1)))
   # par(mfrow=c(1,2))
    image(ker,col=r)
    
    ybin=numeric(0)
    xbin=numeric(0)
    
    for(l in 1:nrow(cuadri1)) {
    ybin[l]=which.min(abs(ker$y-y[l]))-1
    xbin[l]=which.min(abs(ker$x-x[l]))-1
    
    }
              
  #phist=hist(cuadri1$lptn,breaks=seq(0,1,1/25))

  #image(ker, col=r)
  #write.table(ker$z,paste('histogramas2d/hist',formatC(k, width=3, flag="0"),sep='_'),row.names=FALSE,col.names=FALSE)
  write.table(cbind(cuadri1$pt,xbin,ybin),paste('eventosfin/ev',formatC(k, width=3, flag="0"),sep='_'),row.names=FALSE,col.names=FALSE)
  
  # Adjust binning (interpolate - can be computationally intensive for large datasets)
  
  
  #write.table(cuadri1,paste('eventos/event',formatC(k, width=3, flag="0"),sep='_'),row.names=FALSE,col.names=FALSE)
}

d=matrix(0,nrow=nev,ncol=nev)
for(k in 1:nev) {
  if(k==nev) {
    break
  }
  for(p in (k+1):nev) {
    
    
  
    
    ind1=(estables$ev==p)
    ind2=(estables$ev==k)
    
    cuadri1=estables[ind1,c('eta','phi','ips')]
    cuadri2=estables[ind2,c('eta','phi','ips')]
    
    n=nrow(cuadri1)
    m=nrow(cuadri2)
    
  
    castres=as.matrix(parDist(as.matrix(rbind(cuadri1,cuadri2))))
    
    custos=t(castres[-seq(1,n),-seq(n+1,m+n)])
    d[k,p]=max(custos)
    
    
  }
}
max(d)

library(MASS)
# Default call 

r=100
results=matrix(0,ncol=4,nrow=r)

  
indtest1=(estables$ev==sample(1:200,1))
indtest2=(estables$ev==sample(200:400,1))
indtest3=(estables$ev==sample(400:600,1))
indtest4=(estables$ev==sample(600:800,1))




c1=estables[indtest1,c('ips','pt','rapi','phi')]
c2=estables[indtest2,c('ips','pt','rapi','phi')]
c3=estables[indtest3,c('ips','pt','rapi','phi')]
c4=estables[indtest4,c('ips','pt','rapi','phi')]

i03=c3$ips[c3$ips>0]
i04=c4$ips[c4$ips>0]

p03=c3$pt[c3$ips>0]
p04=c4$pt[c4$ips>0]

# i03=c3$ips[c3$ips>0 & c3$ips<10]
# i04=c4$ips[c4$ips>0 & c4$ips<10]
# 
# p03=c3$pt[c3$ips>0 & c3$ips<10]
# p04=c4$pt[c4$ips>0 & c4$ips<10]

sum(c3$ips>0)
sum(c4$ips>0)

plot(i04~log(p04),col='red')
points(i03~log(p03),col='purple')


plot(estables$mag_pv[estables$mag_pv<1e3])



i03=c3$ips[c3$ips>0.5]
i04=c4$ips[c4$ips>0.5]

p03=c3$pt[c3$ips>0.5]
p04=c4$pt[c4$ips>0.5]

cor(i03,p03)
cor(i04,p04)
par(mfrow=c(1,2))
plot(i03~ p03)
abline(lm(i03~ p03))
plot(i04~ p04)
abline(lm(i04~ p04))


res1=kgroups(i03, 2, iter.max = 30, nstart = 5, cluster = NULL)
res2=kgroups(i04, 2, iter.max = 30, nstart = 5, cluster = NULL)


#res$cluster
results[i,1]=mean(i03[res1$cluster==1])
results[i,2]=mean(i03[res1$cluster==2])
results[i,3]=mean(i04[res2$cluster==1])
results[i,4]=mean(i04[res2$cluster==2])


sum((results[,1]+results[,2])>(results[,3]+results[,4]))

par(mfrow=c(2,2))
plot(i03[res1$cluster==1],col='red')
plot(i03[res1$cluster==2])

plot(i04[res2$cluster==1],col='red')
plot(i04[res2$cluster==2])



plot(c3$ips[c3$ips>0])

# tomar so particulas cargadas con production vertex 
# menor que 1m. senón estamos disolvendo o efecto das B
# coller so os ip con pts altos. senon poden ser ruido do weaksingleboson. 

# conseguir un plot no que se vexan moi separados 
# en cada esquina un evento. 






plot(c3$ips,col='blue')
plot(c4$ips)


plot(c1$ips,col='red')
plot(c2$ips,col='green')




plot(density(c1$ips[c1$ips>10]),type='l')
plot(density(c1$ips[c1$ips>10]),type='l')
lines(density(c2$ips[c2$ips>10]),type='l',col='red')
plot(c3$ips)



distancias <- read.table("~/particarlos/pythia8303/rpythia/array.txt", row.names=NULL, quote="\"", comment.char="")

dis1=read.table("~/particarlos/pythia8303/rpythia/disthaler.txt", row.names=NULL, quote="\"", comment.char="")
dis2=read.table("~/particarlos/pythia8303/rpythia/distD.txt", row.names=NULL, quote="\"", comment.char="")

write.table(d,file="distancias.dat",row.names=FALSE,col.names=FALSE)
#distancias <- read.table("distancias800.dat", quote="\"", comment.char="")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
d=dis1
sim=as.matrix(d+t(d))/2
    
diag(sim)=0


  
  perm=sample(1:600,600,replace = FALSE, prob = NULL)
  shuffle=sim[perm,perm]
  
  library(energy)
  
  
  
  #simd=as.dist(sim)
  simp=as.dist(shuffle)
 # simd2=as.dist(sim[(:),(1:473)])
  
  # https://github.com/mariarizzo/k groups do 2019
  res=kgroups(simp, 3, iter.max = 300, nstart = 200, cluster = NULL)
  resp=kgroups(x=as.dist(sim), k=4, iter.max = 1000, nstart = 4000, cluster = NULL)
  tags=resp$cluster
  tags
  tags2=tags
  
  
  #res2=kgroups(simd2, 2,iter.max = 15, nstart = 1, cluster = NULL)
  
  #order(perm)
  tags2=res$cluster[order(perm)]
  tags2
  
  
  step=floor(length(tags2)/3)
  getmode(tags2[1:step])
  table(tags2[1:step])
  getmode(tags2[(step+1):(2*step)])
  table(tags2[(step+1):(2*step)])
  getmode(tags2[(2*step+1):(3*step)])
  table(tags2[(2*step+1):(3*step)])
  getmode(tags2[(3*step+1):(4*step)])
  table(tags2[(3*step+1):(4*step)])
  
  
  step=floor(length(tags)/4)
  getmode(tags[1:step])
  getmode(tags[(step+1):(2*step)])
  getmode(tags[(2*step+1):(3*step)])
  getmode(tags[(3*step+1):(4*step)])
  
  table(tags)


variables=c('le','lmu','ltau','b','q','s')


unique(estables$name)




library(ggplot2)
ggplot(estables, aes(name)) + geom_bar()

  plot(density(estables$pt[(abs(estables$pt)<2)& (estables$name=='mu+' | estables$name=='mu-')]))


  plot(log(ipss[estables$proc=='t']+1))
  plot(log(ipss[estables$proc=='t']+1))
  
  plot(estables$ips[estables$proc==''])
  lines(density(estables$et[estables$proc=='t']))
  plot(density(estables$et[estables$proc=='w']))
  
  
  
  plot(density(log(estables$pt[estables$proc=='t'])))
  
  #prob=(cuadri1$n1)^2-cuadri1$n2^2-cuadri1$n3^2-cuadri1$n4^2
  #probdif=(cuadri1$n1[1]-cuadri2$n1[1])^2-(cuadri1$n2[1]-cuadri2$n2[1])^2-(cuadri1$n3[1]-cuadri2$n3[1])^2-(cuadri1$n4[1]-cuadri2$n4[1])^2




mean(ipss[estables$proc=='t'])/sd(ipss)
mean(ipss[estables$proc=='s'])/sd(ipss)
mean(ipss[estables$proc=='h'])/sd(ipss)
mean(ipss[estables$proc=='w'])/sd(ipss)

mean(estables$n1[estables$proc=='t'])
mean(estables$n2[estables$proc=='s'])
mean(estables$n3[estables$proc=='h'])
mean(estables$n4[estables$proc=='w'])



emedias=aggregate(estables$e, list(estables$ev), mean)
ptmedios=aggregate(estables$pt, list(estables$ev), mean)
ipmedios=aggregate(estables$ips, list(estables$ev), mean)
etamedios=aggregate(estables$eta, list(estables$ev), mean)

plot(etamedios$x~phimedios$x)
points(etamedios$x[lista$x=='w']~phimedios$x[lista$x=='w'],col='red',pch=1) 
points(etamedios$x[lista$x=='h']~phimedios$x[lista$x=='h'],col='green',pch=1) 



ymedios=aggregate(estables$rapi, list(estables$ev), mean)
pv=sqrt(estables$x^2+estables$y^2)
pvmedios=aggregate(pv, list(estables$ev), mean)
phimedios=aggregate(estables$phi, list(estables$ev), mean)

plot(ptmedios)
ind=1:nev
points(ymedios$x[lista$x=='t']~ind[lista$x=='t'],col='brown',pch=20)
points(ymedios$x[lista$x=='w']~ind[lista$x=='w'],col='purple',pch=20)
points(ymedios$x[lista$x=='h']~ind[lista$x=='h'],col='orange',pch=20)
points(ymedios$x[lista$x=='s']~ind[lista$x=='s'],col='cyan',pch=20)

plot(ipmedios)
plot(ptmedios)
plot(ymedios)
ind=1:nev
points(ptmedios$x[lista$x=='t']~ind[lista$x=='t'],col='brown',pch=20)
points(ptmedios$x[lista$x=='w']~ind[lista$x=='w'],col='purple',pch=20)
points(ptmedios$x[lista$x=='h']~ind[lista$x=='h'],col='orange',pch=20)
points(ptmedios$x[lista$x=='s']~ind[lista$x=='s'],col='cyan',pch=20)

plot(ipmedios)
ind=1:nev
points(ipmedios$x[lista$x=='t']~ind[lista$x=='t'],col='brown',pch=20)
points(ipmedios$x[lista$x=='w']~ind[lista$x=='w'],col='purple',pch=20)
points(ipmedios$x[lista$x=='h']~ind[lista$x=='h'],col='orange',pch=20)
points(ptmedios$x[lista$x=='s']~ind[lista$x=='s'],col='cyan',pch=20)





pvmedios=aggregate(pv, list(estables$ev), mean)

plot(log(ptmedios$x[lista$x!='h']),log(emedias$x[lista$x!='h']),pch=20,xlab='log(mean(p_t))',ylab='log(mean(E))')
points(log(ptmedios$x[lista$x=='t']),log(emedias$x[lista$x=='t']),col='red',pch=20)
points(log(ptmedios$x[lista$x=='s']),log(emedias$x[lista$x=='s']),col='cyan',pch=20)
points(log(ptmedios$x[lista$x=='w']),log(emedias$x[lista$x=='w']),col='purple',pch=20)
points(log(ptmedios$x[lista$x=='h']),log(emedias$x[lista$x=='h']),col='darkorange')

plot(log(ptmedios$x),log(emedias$x))
points(log(ptmedios$x[tags2=='1']),log(emedias$x[tags2=='1']),col='red')
points(log(ptmedios$x[tags2=='2']),log(emedias$x[tags2=='2']),col='blue')
points(log(ptmedios$x[tags2=='3']),log(emedias$x[tags2=='3']),col='green')
points(log(ptmedios$x[tags2=='4']),log(emedias$x[tags2=='4']),col='darkorange')



plot(log(pvmedios$x),log(ipmedios$V1))
points(log(pvmedios$x[lista$x=='t']),log(ipmedios$V1[lista$x=='t']),col='red')
points(log(pvmedios$x[lista$x=='w']),log(ipmedios$V1[lista$x=='w']),col='green')
points(log(pvmedios$x[lista$x=='h']),log(ipmedios$V1[lista$x=='h']),col='darkorange')
points(log(pvmedios$x[lista$x=='s']),log(ipmedios$V1[lista$x=='s']),col='blue')

plot((ptmedios$x),(pvmedios$x))
points((ptmedios$x[lista$x=='t']),((pvmedios$x[lista$x=='t'])),col='red')
points((ptmedios$x[lista$x=='w']),(pvmedios$x[lista$x=='w']),col='green')
points(log(ptmedios$x[lista$x=='h']),log(ipmedios$V1[lista$x=='h']),col='darkorange')
points(log(ptmedios$x[lista$x=='s']),log(ipmedios$V1[lista$x=='s']),col='blue')

plot(density(log(ipmedios$V1[lista$x=='t'])),col='red')
lines(density(log(ipmedios$V1[lista$x=='w'])),col='green')
lines(density(log(ipmedios$V1[lista$x=='h'])),col='darkorange')
lines(density(log(ipmedios$V1[lista$x=='s'])),col='blue')




plot(ptmedios$x,pvmedios$x)



library(reticulate)
py_install("Wasserstein")


import('Wasserstein')
use_condaenv("r-reticulate")
conda_create("r-reticulate")

# install SciPy
conda_install("r-reticulate", "Wasserstein")


