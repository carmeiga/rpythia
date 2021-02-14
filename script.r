s=mesmo_proceso('s',0)
h=mesmo_proceso('h',1)
w=mesmo_proceso('w',2)
t=mesmo_proceso('t',3)



estables=rbind(s,h,w,t)
normalizacions=aggregate(estables$e, list(estables$ev), sum)
estables$en=estables$e
nev=nrow(normalizacions)
for(i in 1:nev) {
  estables$en[estables$ev==i]=estables$e[estables$ev==i]/normalizacions$x[i]
}

aggregate(estables$en, list(estables$ev), sum)

d=matrix(0,nrow=nev,ncol=nev)
for(k in 13:nev) {
  for(p in 1:nev) {

proba1=estables[estables$ev==p,c((11:15),22)]
proba2=estables[estables$ev==k,c((11:15),22)]



n=nrow(proba1)
m=nrow(proba2)
custos=matrix(0,nrow=nrow(proba1),ncol=nrow(proba2))

minkowski<-function(cuadri1,cuadri2) {
  return(-cuadri1[4]*cuadri1[4]+sum(cuadri1[1:3]*cuadri1[1:3])-cuadri2[4]*cuadri2[4]+sum(cuadri2[1:3]*cuadri2[1:3])-2*(-cuadri1[4]*cuadri2[4]+sum(cuadri1[1:3]*cuadri2[1:3])))
}

for(i in 1:n) {
  for(j in 1:m){
    
    custos[i,j]=minkowski(as.numeric(proba1[i,]),as.numeric(proba2[j,]))
}
}


d[k,p]=sqrt(wasserstein(proba1$en, proba2$en, p=2, tplan=NULL, costm=abs(custos),prob=TRUE)) # 4-Wasserstein


  }  
}
  

write.table(d,file="distancias.txt",row.names=FALSE,col.names=FALSE) 

distancias <- read.table("~/particarlos/pythia8303/rpythia/distancias.txt", quote="\"", comment.char="")

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









