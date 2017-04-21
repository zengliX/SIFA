#### assisting functions ############

### plot genotype matrix #####
genotype_plot=function(Z,xname='subclone id',yname='mutation id',titlename=NULL)
{
 rownames(Z)=colnames(Z)=c()

 Z=as.matrix(Z)
 temp=melt(Z) 
 p <- ggplot(temp, aes(X2,X1,fill=value) ) + geom_tile() + 
   scale_fill_gradient(low = "white", high = "red")+
   xlab(xname) +
   ylab(yname) +
   ggtitle(titlename) +
   theme(plot.title = element_text(hjust = .5))
 p
}

### plot copy number matrix ######
CN_plot=function(L,xname='subclone id',yname='mutation id',titlename=NULL)
{
  rownames(L)=colnames(L)=c()

  L=as.matrix(L)
  temp=melt(L) 
  p <- ggplot(temp, aes(X2,X1,fill=value) ) + geom_tile() +
        scale_fill_gradient2(low = "green", high = "red",mid='black',midpoint = 2,limit=c(0,4)) +
    theme_minimal() +
  xlab(xname) +
    ylab(yname) +
    ggtitle(titlename) + 
    theme(plot.title = element_text(hjust = .5))
  p
}


### plot prediction errors ######

err_plot=function(E,xname='samples',yname='error',titlename=NULL)
{
  marg=(max(E)-min(E))/10
  range2=c(min(E)-marg,max(E)+marg)
  rownames(E)=colnames(E)=c()
  E=as.matrix(E)
  
  temp=melt(E)
  p= ggplot(temp,aes(X2,value,group=X1)) + geom_line() +ylim(range2) +
    xlab(xname) + ylab(yname) + ggtitle(titlename) + 
    theme(plot.title = element_text(hjust = .5))
  p
}




##### VAF vs total reads ##########
pD_plot=function(X,D,folder='.',name='')  # plot VAF against total reads for each sample
{
  P=X/D
  
  temp_fun=function(a,b)
  {
   plot(a,b,xlab = 'total reads',ylab = 'VAF') 
  }
  
  l=ncol(D)
  pdf(paste(folder,'/',name,'VAF_vs_D.pdf',sep=''))
  for(i in 1:l)
  {
    temp=temp_fun(D[,i],P[,i])
    temp
  }
  dev.off()
}




################################################
################ fraction change plot ########
################################################

frac_plot=function(Frac,title=NULL)
{

Frac2=as.data.frame(cbind(c(1:nrow(Frac)),Frac))

if(!is.null(colnames(Frac))[1]){
  colnames(Frac2)=c('subclone',colnames(Frac))
} else {
  colnames(Frac2)=c('subclone',paste('sample',c(1:ncol(Frac))))
}


Frac2=gather(Frac2,key = 'sample',value = 'fraction',c(2:ncol(Frac2)))
Frac2$subclone=as.factor(Frac2$subclone)

p=ggplot(data=Frac2) + geom_area(aes(x=sample,y=fraction,fill=subclone,group=subclone),
                            position = 'stack')+ggtitle(title) + 
                  theme(plot.title = element_text(hjust = .5))

p
}


################################################
################ grouped line plot ########
################################################
group_line_plot=function(line_mat,group_vec,title=NULL,xlab='sample',ylab='VAF',sep=FALSE)
{
  if(!is.null(colnames(line_mat))[1]){
    samp_names=colnames(line_mat)
  } else {
    samp_names=paste("sample",c(1:ncol(line_mat)),sep='')
  }
  
  colnames(line_mat)=samp_names
  temp=as.data.frame(cbind(group_vec,c(1:nrow(line_mat)),line_mat))
  colnames(temp)[1:2]=c('group','mut')
  temp$group=as.factor(temp$group)
  temp=gather(temp,key='sample',value='VAF',c(3:ncol(temp)))
  
  if(sep){
    l=length(unique(group_vec))
    pic= ggplot(data = temp) + geom_line(aes(x=sample,y=VAF,group=mut,color=group)) + 
      facet_wrap(~group,nrow = ceiling(l/2)) +
     xlab(xlab)+ylab(ylab)+ggtitle(title)
  } else {
  pic=ggplot(data=temp)+geom_line(aes(x=sample,y=VAF,group=mut,color=group))+xlab(xlab)+ylab(ylab)+ggtitle(title)
  }
  pic
}


################################################
################ tree plot ####################
################################################

tree_plot=function(tree)
{
  temp=tree
  temp[1]=1
  tree_df=data.frame(parent=temp,id=c(1:length(temp)))
  tree_graph=graph.data.frame(tree_df)
  plot(tree_graph, layout=layout.reingold.tilford,main="tree structure")
}


###########################################################
############ total reads plot ############################
############################################################

total_r_plot=function(D,log_r=F)
{
  if(!is.null(colnames(D))[1]){
    samp_names=colnames(D)
  } else {
    samp_names=paste("sample",c(1:ncol(D)),sep='')
  }
  
  D2=as.data.frame(D)
  colnames(D2)=samp_names
  
  if(log_r){ # take log of reads
    D2=D2+1
    m=apply(D2,2,mean)
    D2=as.data.frame( log(as.matrix(D2)%*%diag(1/m)))
    colnames(D2)=samp_names
  }
  
  D2$points=c(1:nrow(D))
  
  plot_frame=gather(D2,key="sample",value="reads",1:ncol(D))
  
  l=ncol(D)
  ggplot(plot_frame) + geom_point(aes(x=points,y=reads,color=sample)) + 
    facet_wrap(~sample,nrow = ceiling(l/2))
}





##################################################################
############## get samples for given parameters ##################
##################################################################

get_samp=function(MCMCout,var_name,subset=NULL)  
{
  #subset: a subset of samples
  
  if(is.null(subset)){subset=c(1:length(MCMCout))}
  MCMCout=MCMCout[subset]
  
  N=length(MCMCout)  # number of samples
  temp=MCMCout[[1]][[var_name]]
  if(is.vector(temp)){temp=matrix(temp,nrow = 1)}
  
  d1=nrow(temp)  
  d2=ncol(temp)
  
  if(is.null(d1))#no row dimension,then a scalar
  {
    out=c(1:N)
    for(z in 1:N)
    {
      temp=MCMCout[[z]][[var_name]]
      out[z]=temp
    }
  } else if(d1==1)  # if parameter is vector, get samples as a matrix
  {
    out=matrix(0,N,d2)
    for(z in 1:N)
    {
      temp=MCMCout[[z]][[var_name]]
      out[z,]=temp
    }
  } else {  # else get samples as an array
    out=array(0,dim=c(d1,d2,N))
    for(z in 1:N)
    {
      temp=MCMCout[[z]][[var_name]]
      out[,,z]=temp
    }
  }
  out
}

################################################################
################ evaluate one posteriror sample ################
################################################################

clone_sort=function(Z1,Z0,Params)  # pick permutation with minimum error on Z
{
  #Z1: sample mutation matrix
  #Z0: true mutation matrix
  a=Params$K
  ords=permutations(a-1,a-1,v=2:a) # all possible permutations excluding normal subclone
  ords=cbind(1,ords)
  
  l=nrow(ords) # number of possible orders
  
  loc=0
  v=Inf
  for(i in 1:l)
  {
    Znew=Z1[,ords[i,]]
    temp=sum(abs(Znew-Z0)) # measure of this permutation
    if(temp<v){
      v=temp
      loc=i
    }
  }
  
  v=v/(Params$K*Params$J)
  
  out=list(measure=v,ord=ords[loc,])
  out
}


############################################################################
################ get samples with different tree structures ################
############################################################################

out_trees=function(MCMCout)
{
  l=length(MCMCout)
  all_tree=get_samp(MCMCout,'Ttree')
  
  all_tree=all_tree[!duplicated(all_tree),,drop=F]
  
  out=all_tree
  out
}

## get sample numbers for one tree structure
onetree_samp=function(MCMCout,tree)
{
  l=length(MCMCout)
  
  a=c(1:l)
  for(i in 1:l)
  {
    temp=MCMCout[[i]]$Ttree
    if(all(temp==tree)){a[i]=1} else{a[i]=0}
  }
  
  out=which(a==1)
  out
}


### change tree after reordering subclones
get_newtree=function(old_tree,ord)
{
 K = length(old_tree) 
 out = old_tree
 for(i in 2:K){
   out[i] = which(ord == old_tree[ord[i]])
 }
 out
}


tree_check=function(tree)
{
  out=T

  for(i in 1:length(tree)){
    if(i==1 && tree[1]!=0){return(F)}
    if(i==2 && tree[2]!=1){return(F)}
    if(tree[i]>=i){return(F)}
  }
  out
}

#############################################################
################ DISPLAY PARAMETERS ################
#############################################################
par_disp = function(Params,MCMC_par)
{
  cat("Input information:\n")
  cat("Number of sequencing samples:",Params$T0,"\n")
  cat("Number of loci:", Params$J,"\n")
  cat("Number of genome segments:", nrow(Params$segments),"\n")
  cat("\n")
  
  cat("Sampling parameters:\n")
  cat("Tuning samples:",MCMC_par$Ntune,"\n")
  cat("Burn-in samples:",MCMC_par$burnin,"\n")
  cat("Number of samples to be kept:", MCMC_par$Nsamp,"\n")
  cat("Number of MCMC chains:",MCMC_par$Nchain,
      ", with temperature increment:",MCMC_par$delta_T,"\n")
}



