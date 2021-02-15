s=mesmo_proceso('s',0)
h=mesmo_proceso('h',1)
w=mesmo_proceso('w',2)
t=mesmo_proceso('t',3)

library(transport)

## user-defined distance function
library(RcppArmadillo)
# Use RcppXPtrUtils for simple usage of C++ external pointers
library(RcppXPtrUtils)
library(parallelDist)
## user-defined distance function
minkowskiFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
  return (arma::accu(arma::square(A - B))-2*pow(A(0,0)-B(0,0),2));
}", depends = c("RcppArmadillo"))


# distance matrix for user-defined euclidean distance function (note that method is set to "custom")




estables=rbind(s,h,w,t)
normalizacions=aggregate(estables$e, list(estables$ev), sum)
estables$en=estables$e
nev=nrow(normalizacions)
for(i in 1:nev) {
  estables$en[estables$ev==i]=estables$e[estables$ev==i]/normalizacions$x[i]
}

aggregate(estables$en, list(estables$ev), sum)

d=matrix(0,nrow=nev,ncol=nev)
for(k in 1:nev) {
  for(p in 1:nev) {

proba1=estables[estables$ev==p,c(22,14,(11:13))]
proba2=estables[estables$ev==k,c(22,14,(11:13))]



n=nrow(proba1)
m=nrow(proba2)
custos=matrix(0,nrow=nrow(proba1),ncol=nrow(proba2))


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

res=parDist(as.matrix(rbind(proba1[,-1],proba2[,-1])), method="custom", func = minkowskiFuncPtr)
# }

castres=as.matrix(res)



custos=t(castres[-seq(1,n),-seq(n+1,m+n)])


d[k,p]=sqrt(wasserstein(proba1$en, proba2$en, p=2, tplan=NULL, costm=abs(custos),prob=TRUE)) # 4-Wasserstein


  }  
}
  


write.table(d,file="distancias.dat",row.names=FALSE,col.names=FALSE)
distancias <- read.table("distancias.dat", quote="\"", comment.char="")

sim=as.matrix(distancias+t(distancias))/2

diag(sim)=0

library(energy)

simd=as.dist(sim)

# https://github.com/mariarizzo/kgroups do 2019
res=kgroups(simd, 4, iter.max = 10, nstart = 1, cluster = )



unique(estables$name)




library(ggplot2)
ggplot(estables, aes(name)) + geom_bar()

  plot(density(estables$pt[(abs(estables$pt)<2)& (estables$name=='mu+' | estables$name=='mu-')]))













