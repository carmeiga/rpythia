s=mesmo_proceso('t',0)
h=mesmo_proceso('h',1)
w=mesmo_proceso('w',2)
t=mesmo_proceso('t',3)



estables=rbind(s,h,w,t)
normalizacions=aggregate(estables$e, list(estables$ev), sum)
normalizado=estables

for(i in 1:nrow(normalizacions)) {
  normalizado$e[estables$ev==i]=estables$e[estables$ev==i]/normalizacions$x[i]
}

aggregate(normalizado$e, list(normalizado$ev), sum)

proba1=normalizado[normalizado$ev==1,(11:15)]
proba2=normalizado[normalizado$ev==2,(11:15)]

n=nrow(proba1)
m=nrow(proba2)
custos=matrix(0,nrow=nrow(proba1),ncol=nrow(proba2))

minkowski<-function(cuadri1,cuadri2) {
  return(-cuadri1[5]^2-cuadri2[5]^2-2*(-cuadri1[4]*cuadri2[4]+sum(cuadri1[1:3]*cuadri2[1:3])))
}

for(i in 1:n) {
  for(j in 1:m){
    
    custos[i,j]=minkowski(as.numeric(proba1[i,]),as.numeric(proba2[j,]))
}
}

wasserstein(proba1$e, proba2$e, p=2, tplan=NULL, costm=abs(custos),prob=TRUE)


  
  





unique(estables$name)




library(ggplot2)
ggplot(estables, aes(name)) + geom_bar()

  plot(density(estables$pt[(abs(estables$pt)<2)& (estables$name=='mu+' | estables$name=='mu-')]))









