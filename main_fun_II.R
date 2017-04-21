##################################################
##### parameter random initiation ###############
##################################################

init_gen=function(Params,tree=NULL)
{ # phi: read depth
  
  out=list()
  
  
  #initial tree
  if(is.null(tree)){
    out$Ttree=c(0,1,rep(1,Params$K-2))
  } else {
    out$Ttree=tree
  }
  
  if(Params$K>2){
    for(z in 3:Params$K){ out$Ttree[z]=sample(1:(z-1),1) }
  }
  
  # initial pi
  out$pi=0.5
  
  # initial L
  tempL=tempZ=matrix(1,Params$J,Params$K)
  segs=Params$segments
  for(i in 1:nrow(segs))
  {
    cur_seg=seq(segs[i,1],segs[i,2])
    temp=propose_L_cpp(out$Ttree,tempZ[cur_seg,,drop=F],Params)
    for(j in cur_seg){
    tempL[j,]=temp
    }
  }
  out$L=tempL
  
  # initial Z
  tempZ=matrix(0,Params$J,Params$K)
  for(j in 1:Params$J)
  {
    tempZ[j,]=propose_Z_cpp(out$Ttree,tempL[j,],Params)
  }
  out$Z=tempZ
  
  # initial theta
  out$theta=matrix(rgamma(Params$K*Params$T0,2,1),Params$K,Params$T0)
  
  # initial F
  temp=apply(out$theta,2,sum)
  out$Frac=out$theta%*%diag(1/temp)
  
  # initial M
  out$M=out$L %*% out$Frac
  
  out
}


##################################################
##### generate pars from prior distr ###############
##################################################

prior_gen=function(Params)
{ 
  out=list()
  
  
  #initial tree
    out$Ttree=c(0,1,rep(1,Params$K-2))
  if(Params$K>2){
    for(z in 3:Params$K){ out$Ttree[z]=sample(1:(z-1),1) }
  }
  
  # initial pi
  out$pi= rbeta(1,Params$a_pi,Params$b_pi)
  
  # initial Z
  tempL=tempZ=matrix(1,Params$J,Params$K)
  
  tempZ=matrix(0,Params$J,Params$K)
  w = Params$zeta^c(1:Params$max_mut)
  w= w/sum(w)
  for(j in 1:Params$J)
  {
    if(Params$K > 2){
      temp = c( sample(c(2:Params$K),1) ,sample(c(1:Params$max_mut),1,prob = w) )
    } else {
      temp = c(2, sample(c(1:Params$max_mut),1,prob = w) )
    }
    
    tempZ[j,]= Zo_to_Z_cpp(t(as.matrix(temp)),out$Ttree)
  }
  out$Z=tempZ
  
  # initial L
  segs=Params$segments
  for(i in 1:nrow(segs))
  {
    cur_seg=seq(segs[i,1],segs[i,2])
    without_cnv = rbinom(1,1,out$pi)
    if(without_cnv == 1) {
      temp = rep(2,Params$K)
    } else {
      a = max(out$Z[cur_seg,])
      b = c((a-2):(Params$max_CN -2) )
      b = setdiff(b,0)
      if(Params$K > 2){
        temp = c(sample(c(2:Params$K),1),sample(b,1))
      } else {
        temp = c(2,sample(b,1))
      }
      temp = Lo_to_L_cpp(t(as.matrix(temp)),out$Ttree)
    }
    
    for(j in cur_seg){
      tempL[j,]=temp
    }
  }
  out$L=tempL
  
  # initial theta
  out$theta=matrix(rgamma(Params$K*Params$T0,Params$r,1),Params$K,Params$T0)
  
  # initial F
  temp=apply(out$theta,2,sum)
  out$Frac=out$theta%*%diag(1/temp)
  
  # initial M
  out$M=out$L %*% out$Frac
  
  out
}


#######################################################
########## heuristic init gen ####################
#######################################################
heuristic_tree_gen=function(X,D,K)
{

 P=X/D # observed VAF
 P_cluster=kmeans(P,K-1)
 mut_g=P_cluster$cluster
 centers=P_cluster$centers
 
 # initialize tree
 tree=c(0:(K-1))
 vaf_mean=apply(centers,1,mean)
 ord=order(vaf_mean,decreasing = T)
 centers=centers[ord,]
 centers=rbind(1,centers)
 
 for (i in 3:K){
   a=c()
   for (j in 1:(i-1)){
      a[j]=sum( centers[j,] >= centers[i,]  )
   }
   tree[i]=sample(1:I(i-1),1,prob = a)
 }
 
 tree
}



#######################################################
########## acquire one MCMC sample ####################
#######################################################

MCMC_onesamp=function(X,D,Params,init,adap,temper=1)
{ # init: last round MCMC sample
  # temper: temperature
  # adap: adaptive tuning parameters
  
  out=init
  
  
  ########## update pi #####
  
  out$pi=samp_pi(out,Params,temper)
  
  ########## update F (or Theta)
  
  out$theta=samp_theta_all(out,Params,X,D,adap$theta_tune,temper)
  temp=apply(out$theta,2,sum)
  out$Frac=out$theta%*%diag(1/temp) # update F with new theta values
  
  
  ######### update L matrix ##########
  out$L=samp_L_all(out,X,D,Params,temper)
  
  ###### update M matrix #########
  out$M=out$L %*% out$Frac
  
  #### update Z matrix #########
  out$Z=samp_Z_all(out,X,D,Params,temper)
  
  ##### update tree structure ########
  
  if(Params$K>2)
  {
    Lo=L_to_Lo_cpp(out$L)
    Zo=Z_to_Zo_cpp(out$Z)
    
    u=runif(1)
    if(u<0.15){ 
      slice_out=slice_samp_tree(X,D,out,Lo,Zo,Params,temper=temper)
      if(!is.null(slice_out))
      {
        out$Ttree=slice_out$Ttree
        out$L=slice_out$L
        out$Z=slice_out$Z
      }
    } else {
      out=MH_samp_tree(out$Ttree,X,D,out,Params,Lo,Zo,temper)
    }
  }
  
  ###### keep track of likelihood ########
  p=prob_XD(X,D,out$Frac,out$L,out$Z,Params$phi)
  out$likelihood=p
  
  out
}

























