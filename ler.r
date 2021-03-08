mesmo_proceso<-function(pro,tot) {
  setwd("~/particarlos/pythia8303/rpythia")
  pat=paste("_",pro,".txt$",sep="")
  
  list_of_files <- list.files(path = ".", recursive = TRUE,
                              pattern = pat, 
                              full.names = FALSE)
  
  

  separ=as.numeric(read.table(list_of_files[5], quote="\"", comment.char="")$V1)
  
  n=length(separ)
  ev=n/2-1
  
  
  library(readr)
  
  
  
  for(i in 1:ev) {
    
    temp <- read_table2(list_of_files[4], col_names=FALSE,skip =(separ[2*i+1]+2),n_max=(separ[2*i+2]-separ[2*i+1]-5))
    
    temp$ev=rep(i+tot*ev,nrow(temp))
    if (i==1)
      out=temp
    else 
      out=rbind(out,temp)
  }
  
  nomes= c('no','id','name','status', 'mother1','mother2','daughter1','daughter2','colour1','colour2','p_x','p_y','p_z','e','m','ev')
  
  colnames(out)=nomes
  
  pt <- as.numeric(read_table2(list_of_files[3], col_names = FALSE)[1,])
  out$pt=pt[-length(pt)]
  
  x <- as.numeric(read_table2(list_of_files[7], col_names = FALSE)[1,])
  out$x=x[-length(x)]
  
  y <- as.numeric(read_table2(list_of_files[8], col_names = FALSE)[1,])
  out$y=y[-length(y)]
  
  z <- as.numeric(read_table2(list_of_files[9], col_names = FALSE)[1,])
  out$z=z[-length(z)]
  
  # q <- as.numeric(read_table2(list_of_files[1], col_names = FALSE)[1,])
  # out$q=q[-length(q)]
  # 
  # spin <- as.numeric(read_table2(list_of_files[6], col_names = FALSE)[1,])
  # out$spin=spin[-length(spin)]
  # 
  et <- as.numeric(read_table2(list_of_files[2], col_names = FALSE)[1,])
  out$et=et[-length(et)]
  # 
  estables_n=c('e-','e+','mu+','mu-','K+','K-','pi+','pi-','p+','pbar-')#,'n0','nbar0','gamma','K_L0') #sen neutras!
  # parenteses indican particulas intermediarias (desintegranse)
  
  estables=subset(out, (name %in% estables_n))
  
  final=cbind(estables,rep(pro,nrow(estables)))
  names(final)[22] <- "proc"
  final$b=as.numeric(final$name=='p+' | final$name=='pbar-' | final$name=='n0' | final$name=='nbar0')
  final$le=as.numeric(final$name=='e+' | final$name=='e-')
  final$lm=as.numeric(final$name=='mu+' | final$name=='mu-')
  
  
  
  
  return(final)
}

