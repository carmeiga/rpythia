setwd("~/particarlos/pythia8303/rpythia")
source('ler.r')


s=mesmo_proceso('s',0)
h=mesmo_proceso('h',1)
w=mesmo_proceso('w',2)
t=mesmo_proceso('t',3)



library(transport)
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





cppFunction('double ip(arma::vec p, arma::vec v) {
  double impact = arma::norm(arma::cross(p,v));
  return impact;
}',depends = c("RcppArmadillo"))

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
normalizacions=aggregate(estables$e, list(estables$ev), sum)
estables$en=estables$e
nev=nrow(normalizacions)
for(i in 1:nev) {
  estables$en[estables$ev==i]=estables$e[estables$ev==i]/normalizacions$x[i]
}

aggregate(estables$proc, list(estables$ev),max)

ipss=ips(cbind(estables$p_x,estables$p_y,rep(0,length(estables$p_x))),cbind(estables$x,estables$y,rep(0,length(estables$p_x))))

estables$ips=ipss

estables$n1=estables$e/estables$et
estables$n2=estables$p_x/estables$et
estables$n3=estables$p_y/estables$et
estables$n4=estables$p_z/estables$et

grupo=as.matrix(cbind(estables$n1,estables$n2,estables$n3,estables$n4))

sigma=cov(grupo)
sigmainv=solve(sigma)

sigmam=chol(sigmainv)
grupo=grupo%*%t(sigmam)


estables$n1=grupo[,1]
estables$n2=grupo[,2]
estables$n3=grupo[,3]
estables$n4=grupo[,4]
normalize=function(x) {
  return((x-min(x))/(max(x)-min(x)))
}



d=matrix(0,nrow=nev,ncol=nev)
for(k in 1:nev) {
  if(k==nev) {
    break
  }
  for(p in (k+1):nev) {
    
   
    ind1=(estables$ev==p)
    ind2=(estables$ev==k)
    
    cuadri1=estables[ind1,c('n1','n2','n3','n4')]
    cuadri2=estables[ind2,c('n1','n2','n3','n4')]
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
#res2=parDist(as.matrix(rbind(disc1,disc2)), method="custom", func = trivial_dist)
if((sum(ip1)<1e-20) | (sum(ip2)<1e-20) ){res3=matrix(rep(0,(n+m)^2),ncol=n+m)
}else{res3=parDist(as.matrix(c(ip1,ip2)),method='mahalanobis')}

castres=sqrt(sqrt(abs(as.matrix(res1)))+(as.matrix(res3)))#+sqrt(as.matrix(res2))



# }





custos=t(castres[-seq(1,n),-seq(n+1,m+n)])


d[k,p]=wasserstein(as.numeric(en1), as.numeric(en2), p=2, tplan=NULL, costm=abs(custos),prob=TRUE)


  }  
}


  


write.table(d,file="distancias.dat",row.names=FALSE,col.names=FALSE)
distancias <- read.table("distancias.dat", quote="\"", comment.char="")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

sim=as.matrix(d+t(d))/2
    
diag(sim)=0


  
  perm=sample(1:nev, nev,replace = FALSE, prob = NULL)
  shuffle=sim[perm,perm]
  
  library(energy)
  
  
  
  #simd=as.dist(sim)
  simp=as.dist(shuffle)
  #simd2=as.dist(sim[(1:473),(1:473)])
  
  # https://github.com/mariarizzo/kgroups do 2019
  #res=kgroups(simd, 4, iter.max = 15, nstart = 1, cluster = NULL)
  resp=kgroups(as.dist(sim), 4, iter.max = 1, nstart = 1, cluster = NULL)
  resp$cluster
  
  
  #res2=kgroups(simd2, 2,iter.max = 15, nstart = 1, cluster = NULL)
  
  #order(perm)
  tags=resp$cluster[order(perm)]
  
  
  step=floor(length(tags)/4)
  getmode(tags[1:step])
  getmode(tags[(step+1):(2*step)])
  getmode(tags[(2*step+1):(3*step)])
  getmode(tags[(3*step+1):(4*step)])


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










