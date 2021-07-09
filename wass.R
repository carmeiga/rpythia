# os primeiros 600 son simulacions. os ultimos 800 son experimentais (sinal+background)
dis1=sqrt(read.table("array.txt", row.names=NULL, quote="\"", comment.char=""))
#
library('plot.matrix')
nobsm=c(1:600,601:750,801:950,1001:1150) #to not introduce bsm higgs
bsm=c(1:600,601:750,801:950,1001:1150,1201:1400) # 1201 to 1250 entries are bsm higgses
#these index have been created to choose not so many weak and top events

dis1=dis1[bsm,bsm]

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


library(FastKNN)

mb=600 # simulated background size
n=650 # real data size
pii=n/(mb+n)

#non hai que separar a mostra!!!!!!!!!
#trax=sample(1:mb,mb) #indices background para adestrar
#traw=sample((mb+1):(nr+mb),n) #indices reais para adestrar

#adestr=dis1[c(trax,traw),c(trax,traw)]
tagadestr=c(rep(0,mb),rep(1,n)) # etiquetas, forman parte da mostra
#tagtest=c(rep(0,1200),rep(1:200)) non os sabemos!!!!

ntot=nrow(dis1)

k=500 # grid of k's on which kfold cv will be used
kf=4 # numero de folds
errtemp=numeric(kf)
require(caret)
flds <- createFolds(tagadestr, k = kf, list = TRUE, returnTrain = FALSE)
names(flds)[1] <- "train"
  
#tagtestsh=tagtest[indsh]
erro1=numeric(k)
  
  #erro2=numeric(k)
  #erro2=0
  
for(r in 1:k) { # r runs over k's grid
   
    
    for(v in 1:kf) {
      indfold=flds[[v]] # indices para validar
      aux=1:ntot
      trainfold <- aux[is.na(pmatch(aux,indfold))] #indices para adestrar
      
      #separation of response
      tagsh=tagadestr[trainfold]
      tagtest=tagadestr[indfold]
      
      ntrain=length(trainfold)
      ntest=length(indfold)
      
      #compute r nearest neighbours for each observation in the fold
      #we cannot use the whole distance matrix! we would have neighbours in the test subset
      
      nn = matrix(0,ntest,r) 
      
      for (i in 1:ntest) {
        disaux=as.matrix(dis1[c(trainfold,indfold[i]),c(trainfold,indfold[i])])
        nn[i,] = order(as.numeric(disaux[(ntrain+1),-(ntrain+1)]))[1:r]
        
      }
      
      tagfinal=matrix(0,nrow=ntest,ncol=r)
      #tagfinal=numeric(n)
      pred=numeric(ntest)
      for(i in 1:ntest) {
        tagfinal[i,]=tagsh[nn[i,1:r]]
        pred[i]=mean(tagfinal[i,])
      }
      nfold=sum(tagtest==1)
      mfold=sum(tagtest==0)
      piifold=nfold/(nfold+mfold)
      hw=pred[tagtest==1]
      hx=pred[tagtest==0]
      errtemp[v]=0.5*((1/length(hx))*sum(hx>piifold)+(1/length(hw))*sum(hw<piifold))
      
    }
   
    erro1[r]=mean(errtemp) 
    
    # erro1[r]=(sum(abs(tagsh-pred)))/1400
    # erro2[r]=sum(abs(pred[1:300]))/300
  }


#plot(erro2,col='red',ylim=c(0,1))
#points(erro1)

#plot(erro2[-1],col='red')
grellaks=(1:length(erro1))
bw0 <- npreg(bws=20,tydat = erro1,txdat=grellaks )
plot(erro1~grellaks,ylim=c(0.43,0.5))
lines(bw0$mean~bw0$eval$grellaks,col='red',lwd=2)

pdf(file = "kfold.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches

plot(erro1~grellaks,ylim=c(0.43,0.5),xlab='k',ylab='Proporci�n de erro')
lines(bw0$mean~bw0$eval$grellaks,col='red',lwd=2)

which.min(bw0$mean)

dev.off()



res=plot(bw0,plot.behavior='plot-data',neval=k)

modparab=lm(erro1~I(grellaks^2)+grellaks)

which.min(res$r1$mean)

plot(erro1,col='black',)
lines(grellaks,res$r1$mean,col='red',lwd=2)

indkopt=which.min(erro1)



# ata aquí escoller o numero de veciños optimo
indkopt=265


theta=numeric(B)
That=numeric(B)
M=numeric(B)
nn = matrix(0,ntot,indkopt)

for (i in 1:ntot)
  nn[i,] = k.nearest.neighbors(i, dis1, k = indkopt)

for(b in 1:1000) {
  outputclas=numeric(ntot)

  tagfinal=matrix(0,nrow=ntot,ncol=indkopt)
  tagrandperm=tagadestr[sample(1:ntot,ntot,replace=FALSE)]
  
  for(i in 1:ntot) {
    tagfinal[i,]=tagrandperm[nn[i,1:indkopt]]
    outputclas[i]=mean(tagfinal[i,])
  }
  
  
  hw=0
  hx=0
  
  hw=outputclas[tagrandperm==1]
  hx=outputclas[tagrandperm==0]
  
  That[b]=2*n*(log((1-pii)/pii)+(1/n)*sum(log(hw/(1-hw))))
  
  indicador=numeric(mb)
  for(i in 1:mb)
    indicador[i]=sum(hw>hx[i])
  
  #indicador=numeric(n)
  # for(i in 1:n)
  #   indicador[i]=sum(hw[i]>hx)
  # sum(indicador)
  
  theta[b]=(1/(n*mb))*sum(indicador)
  M[b]=0.5*((1/mb)*sum(hx>pii)+(1/n)*sum(hw<pii))
}

# os estatísticos sobre a nosa mostra
outputclas=numeric(ntot)
tagfinal=matrix(0,nrow=ntot,ncol=indkopt)

for(i in 1:ntot) {
  tagfinal[i,]=tagadestr[nn[i,1:indkopt]]
  outputclas[i]=mean(tagfinal[i,])
}

hw=0
hx=0
hw=outputclas[tagadestr==1]
hx=outputclas[tagadestr==0]


Thatest=2*n*(log((1-pii)/pii)+(1/n)*sum(log(hw/(1-hw))))

indicador=numeric(mb)
indicador=0
for(i in 1:mb)
  indicador[i]=sum(hw>hx[i])

#indicador=numeric(n)
# for(i in 1:n)
#   indicador[i]=sum(hw[i]>hx)
# sum(indicador)

thetaest=(1/(n*mb))*sum(indicador)
Mest=0.5*((1/mb)*sum(hx>pii)+(1/n)*sum(hw<pii))

par(mfrow=c(1,3))
hist(That)
abline(v=Thatest,col="red")
hist(M,xlim=c(0.40,0.60))
abline(v=Mest,col="red")

hist(theta,xlim=c(0.40,0.60))
abline(v=thetaest,col="red")



thatbsm=Thatest
mestbsm=Mest
thetabsm=thetaest

plot(outputclas[order(indsh)])

# In conclusion, we see that slow permutation method when using the AUC and MCE statis-tics out-performs the other semi-supervised methods and additionally gives comparable per-formance to the supervised methods in detecting the signal in the experimental data



step=floor(length(tags2)/4)
getmode(tags2[1:step])
table(tags2[1:step])
getmode(tags2[(step+1):(2*step)])
table(tags2[(step+1):(2*step)])
getmode(tags2[(2*step+1):(3*step)])
table(tags2[(2*step+1):(3*step)])
getmode(tags2[(3*step+1):(4*step)])
table(tags2[(3*step+1):(4*step)])