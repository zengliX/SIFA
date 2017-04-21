
##################################################################
################ sampling of  tree ################################
##################################################################



#### for slice sampling of tree structure
slice_samp_tree=function(X_f,D_f,samp,Lo,Zo,Params,temper=1,Nmax=70)
{
  # Nmax: max number of samples
  p0=prob_XD(X_f,D_f,samp$Frac,samp$L,samp$Z,Params$phi)
  p0=p0/temper
  logu=p0-rexp(1)
  
  l=length(samp$Ttree)
  for(ct in 1:Nmax)
  {
    new_T= rep(0,Params$K)
    new_T[2]=1
    for(z in 3:l){new_T[z]=sample(c(1:(z-1)),1) } # random generate proposal tree
    
    Lnew=Lo_to_L_cpp(Lo,new_T)
    Znew=Zo_to_Z_cpp(Zo,new_T)
    
    if(any((Lnew-Znew)<0)){next} # make sure proposal is valid
    
    p1=prob_XD(X_f,D_f,samp$Frac,Lnew,Znew,Params$phi)
    p1=p1/temper
    
    if(p1>logu){break}
  }
  
  if(ct==Nmax)
  {out=NULL} else{
    out=list()
    out$Ttree=new_T
    out$L=Lnew
    out$Z=Znew
  }
  
  out
}



##################################################################
##### Metropolis- Hastings sampling ##############################
#####################################################################
MH_samp_tree=function(tree,X_f,D_f,samp,Params,Lo,Zo,temper=1)
{
  p0=prob_XD(X_f,D_f,samp$Frac,samp$L,samp$Z,Params$phi)
  p0= p0 + log_prior_all(samp,Params,temper)

  ###### propose new tree
  l=length(tree)
  
  if(l==3){temp=3} else{ # possible subclone to be rerouted
    temp=sample(c(3:l)) 
  }
  
  for(k in temp) #pick a leaf node to reroute
  {
    if(sum(find_desc_cpp(k,tree))==1){
      child=k;break;
    }
  }
  
  par_set=setdiff(c(1:Params$K),child) # possible set of parent
  par_set=sample(par_set)
  
  for(i in par_set){ # pick a parent
    if(can_append(child,i,tree,samp,Lo,Zo)) {par=i;break;}
  }
  
  
  ######### append child to par
  if(par==tree[child]){ # if appended to the same parent
    return(samp)
  }
  
  # if appended to other node
  new_tree=tree
  new_tree[child]=par
  
  new_theta=samp$theta
  theta_child=samp$theta[child,]
  theta_par=samp$theta[par,]
  theta_min=apply(rbind(theta_child,theta_par),2,min)
  new_theta[child,]=theta_min
  new_theta[tree[child],]=samp$theta[tree[child],]+ theta_child
  new_theta[par,] = samp$theta[par,] - theta_min
  
  if(any(new_theta==0)) {new_theta[new_theta==0]= 10^-4} # add noise to zero elements}
  theta_sum=apply(new_theta,2,sum)
  
  new_frac=new_theta%*%diag(1/theta_sum)
  
  
  ##### sort new tree if violate order constraint
  temp=tree_org(new_tree)
  Lo2=Lo;Zo2=Zo
  new_tree=temp$tree
  
  if(!all(temp$ord==c(1:Params$K)))
  {
    new_frac=new_frac[temp$ord,]
    new_theta=new_theta[temp$ord,]
    
    Lo2[,1]=sapply(Lo[,1],function(x){
      ifelse(x==0,0,which(temp$ord==x))})
    Zo2[,1]=sapply(Zo[,1],function(x){which(temp$ord==x)})
  }
  
  new_samp=samp
  new_samp$Ttree=new_tree
  new_samp$L=Lo_to_L_cpp(Lo2,new_tree)
  new_samp$Z=Zo_to_Z_cpp(Zo2,new_tree)
  new_samp$theta=new_theta
  new_samp$Frac=new_frac
  
  new_samp$M= new_samp$L %*% new_samp$Frac
  new_samp$likelihood= prob_XD(X_f,D_f,new_samp$Frac,new_samp$L,new_samp$Z,Params$phi)
  p1=new_samp$likelihood
  p1= p1 + log_prior_all(new_samp,Params,temper)
    
  ##### calculate acceptance rate 
  acc_prob=min(1,exp(p1-p0))
  if(runif(1)<acc_prob){
    return(new_samp)
  } else {
    return(samp)
  }
}



##################################################################
######### check whether can append #####################
############################################################

can_append=function(child,par,tree,samp,Lo,Zo)
{
  K=length(tree)
  J=nrow(samp$L)
  
  par_L=samp$L[,par]
  par_Z=samp$Z[,par]
  
  for(i in 1:J){
    if(Zo[i,1]==child && Zo[i,2]>par_L[i]) {return(FALSE)}
    if(Lo[i,1]==child && 2+Lo[i,2]<par_Z[i]){return(FALSE)}
  }
  return(TRUE)
}



###############################################################
############## reorganize tree ##############################
############################################################
tree_org=function(tree)
{
  l=length(tree)
  ord=c(1:l)


  ntree=tree
  i=2
  while(i < l)
  {
    parent=ntree[i]
    if(parent>(i-1))
    {
      ord[c(parent,i)]=ord[c(i,parent)]  #keep track of exchanging
      ntree[c(parent,i)]=ntree[c(i,parent)]
      temp=ntree
      ntree[which(temp==parent)]=i
      ntree[which(temp==i)]=parent
      # switch the index of i and its parent
      i=2
    } else {i=i+1}
  }

  out=list()
  out$tree=ntree
  out$ord=ord
  out
}

