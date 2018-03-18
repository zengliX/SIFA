####################################################################
############ get point estimation  ####################
####################################################################

get_point_estimate=function(foldername,rdata)
{
  ### load rdata in this iteration
  cur_rdata=paste(foldername,rdata,sep='/')
  load(cur_rdata)
  # reorganize samples
  chain1 = sample_organize(Trace)
  # handling multiple trees
  all_tree=out_trees(chain1)
  if(nrow(all_tree)>1){
    locs=apply(all_tree,1,function(x){onetree_samp(chain1,x)})  # samples for each tree
    if(length(locs)>3) # only keep top 3 most frequent trees
    { temp=order(sapply(locs,length),decreasing = T)[1:3]
    locs=locs[temp]
    } 
  } else{locs=list(c(1:MCMC_par$Nsamp))}
    
  # for one tree, plot median estimation
  out = list()
  for(j in 1:length(locs))
  {
    point_est = list()
    #extract samples for one tree
    chain2=chain1[locs[[j]]] 
    cur_tree=chain2[[1]]$Ttree
    point_est$Ttree = cur_tree

    # mean F
    temp=get_samp(chain2,'Frac')
    med_F=apply(temp,c(1,2),mean)
    point_est$Frac = med_F
    
    # median L
    temp=get_samp(chain2,'L')
    med_L=apply(temp,c(1,2),median)
    point_est$L = med_L
    
    # median Z
    temp_Z=get_samp(chain2,'Z')
    med_Z=apply(temp_Z,c(1,2),median)
    point_est$Z = med_Z
    
    out[[j]] = point_est
  }
  return(out)
}



