
setwd("~/particarlos/pythia8303/rpythia")
separ=as.numeric(read.table("separadores.txt", quote="\"", comment.char="")$V1)

n=length(separ)
ev=n/2-1


library(readr)



for(i in 1:ev) {
  
  temp <- read_table2("saida.txt", col_names=FALSE,skip =(separ[2*i+1]+2),n_max=(separ[2*i+2]-separ[2*i+1]-5))
 
  temp$ev=rep(i,nrow(temp))
  if (i==1)
    out=temp
  else 
    out=rbind(out,temp)
}

nomes= c('no','id','name','status', 'mother1','mother2','daughter1','daughter2','colour1','colour2','p_x','p_y','p_z','e','m','ev')

colnames(out)=nomes

pt <- as.numeric(read_table2("pt.txt", col_names = FALSE)[1,])
out$pt=pt[-length(pt)]

x <- as.numeric(read_table2("xprod.txt", col_names = FALSE)[1,])
out$x=x[-length(x)]

y <- as.numeric(read_table2("yprod.txt", col_names = FALSE)[1,])
out$y=y[-length(y)]

z <- as.numeric(read_table2("zprod.txt", col_names = FALSE)[1,])
out$z=z[-length(z)]

estables_n=c('e-','e+','mu+','mu-','K+','K-','pi+','pi-','p+','pbar-','n0','nbar0','gamma','K_L0')
# parenteses indican particulas intermediarias (desintegranse)

estables=subset(out, (name %in% estables_n))

unique(estables$name)


library(ggplot2)
ggplot(estables, aes(name)) + geom_bar()

  plot(density(estables$p_x[abs(estables$pt)<1]))









