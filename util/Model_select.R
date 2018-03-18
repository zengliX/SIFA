####################################################################################
################ Bayesian Free Energy ############################################################
################################################################################


BFE_calc=function(X,D,Temperature,onefolder)
{
  #### model selection
  rdata=grep("K.*.Rdata",list.files(onefolder),value = T)
  BayEnergy=list(K=c(),BayEnergy=c())
  
  
  for(i in 1:length(rdata)){
    onerdata=rdata[i]
    temp=paste(onefolder,"/",onerdata,sep = "")
    load(temp)
    
    # calculate likelihood for beta = 0 chain
    p_beta0 = c(1:MCMC_par$Nsamp)
    for(j in 1:MCMC_par$Nsamp){
      temp = prior_gen(Params)
      p_beta0[j] = prob_XD(X,D,temp$Frac,temp$L,temp$Z,Params$phi)
    }
    beta_log_like = rbind(beta_log_like,p_beta0)
    
    
     # free energy
     BayEnergy$K = c(BayEnergy$K,Params$K)
    b_energy = 0
    rev_t = c(1/Temperature,0)
    for(j in 2:length(rev_t)){
      beta_l = beta_log_like[j,]*(rev_t[j-1] - rev_t[j])
      M = max(beta_l)
      beta_l = beta_l - M
      b_energy = b_energy - (log(mean(exp(beta_l))) + M)
    }
    BayEnergy$BayEnergy = c(BayEnergy$BayEnergy, b_energy)
    
  }
  
  plot(BayEnergy$BayEnergy ~ BayEnergy$K,main="Bayes energy selection" ,xlab="K",ylab = "Bayes free enerty")
  save(BayEnergy,file=paste(onefolder,"/BayEnergy_model_selection.Rdata",sep=""))
}

