#####################################################
############ main functon ###########################
######################################################

CallSeg=function(D,mut_loc,foldername)
{
  temp = which(D==0)
  if(length(temp) >0){
    D[temp]=1
  }
  
  cat("calculating CNV segments \n")

  J=nrow(D)
  T0=ncol(D)


  m=apply(D,2,median)
  logD=log(D%*%diag(1/m))
  logD=scale(logD)

  #### create matrix
  Segcall_dat=as.data.frame(cbind(mut_loc[,c('chromosome','position')],logD))
  colnames(Segcall_dat)[3:(2+T0)]=paste("sample",c(1:T0))

  #### run segmentation algorithm
  gamma_range=exp(seq(-1,log(50),length.out = 30))
  BIC=c()
  N=J*T0 # total sample size

  for( cur_gamma in gamma_range){
    cat("calculating BIC for gamma:",cur_gamma,"\n")
    multi_seg=multipcf(data=Segcall_dat,verbose = F, pos.unit="mbp",assembly = "hg19",
                       gamma=cur_gamma,normalize = F)

  
    # calculate BIC
    segs=get_segs(mut_loc, D, multi_seg)
    fitted=get_fit(segs)

    errs=logD-fitted  # fitting errors
    Nsegs=nrow(segs$segs)

    cur_BIC=0
    for(i in 1:T0){
      sigma2=sum(errs[,i]^2)/J
      cur_BIC= cur_BIC + J*log(sigma2)
    }
    cur_BIC=cur_BIC + log(N) * ((T0+1)*Nsegs)

    BIC=c(BIC,cur_BIC)
  }

  
  #### select the gamma with smallest BIC
  opt_gamma=gamma_range[which.min(BIC)]
  multi_seg=multipcf(data=Segcall_dat,verbose = FALSE, normalize = F, gamma = opt_gamma,
                     pos.unit = "mbp")
  segs=get_segs(mut_loc, D, multi_seg)
  fitted=get_fit(segs)


  filename=paste('./',foldername,'/segments.pdf',sep='')
  pdf(filename)
  plot(BIC~gamma_range)
  p=fit_view(logD,fitted)
  plot(p)
  dev.off()


  # return segments
  segs$segs
}


# ############################################################
# ############### assistive function ########################
# ############################################################

temp_fun=function(c,chrom,chr_range)
{
  if(c[1]==chrom && c[2] >= chr_range[1] && c[2] <= chr_range[2] ){
    out=TRUE
  } else{
    out=FALSE
  }
  out
}


###### get segments from multipcf function
get_segs=function(mut_loc, D, multi_seg)
{
  out=list()

  T0=ncol(D)
  J=nrow(D)

  #get segs
  segs=matrix(0,nrow(multi_seg),2)
  for( i in 1:nrow(multi_seg)){
    chr=multi_seg[i,'chrom']
    chr_range=multi_seg[i,c(3:4)]
    loc=which(apply(mut_loc[,1:2],1,temp_fun,chrom=chr,chr_range=chr_range))
    segs[i,]=range(loc)
  }

  val=multi_seg[,c(6:(5+T0))]

  out$segs=segs
  out$val=val

  out
}


##############  get fitted result from get_segs
get_fit=function(segs)
{
  segments=segs$segs
  vals=segs$val

  T0=ncol(vals) #number of samples

  fitted=matrix(0,max(segments),T0)


  for ( j in 1:nrow(segments)){
    temp=seq(segments[j,1],segments[j,2])
    for( i in 1:T0){
      fitted[temp,i]=vals[j,i]
    }
  }

  fitted
}


################ visualize fitted result
fit_view=function(logD,fitted)
{
  J=nrow(logD)
  T0=ncol(logD)
  
  obs=as.data.frame(logD)
  colnames(obs)=paste('sample',c(1:T0))
  obs$group='obs'
  obs$ind=c(1:J)
  
  fitted=as.data.frame(fitted)
  colnames(fitted)=paste('sample',c(1:T0))
  fitted$group='fit'
  fitted$ind=c(1:J)
  
  for_plot=as.data.frame(rbind(obs,fitted))
  for_plot=gather(data = for_plot, key = 'sample',value = 'value',1:T0)
  p= ggplot(data = for_plot) + geom_point(aes(x=ind,y=value,color=group)) +
    facet_wrap(~sample)  + ggtitle("Log total reads fitted by segments")
  
  p
}


