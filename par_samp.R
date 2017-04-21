####### functions for my method ##########


##################################################
#######  Gibbs sampling of phi ##################
##################################################
samp_phi=function(D_f,samp,Params)  # _f stands for variable name in functions
{
  ahat=Params$a_phi + apply(D_f,2,sum)
  bhat=Params$b_phi + apply(samp$M,2,sum)/2

  samps=matrix(rgamma(Params$T0,ahat,bhat),ncol=Params$T0)  # sample from updated parameter

  out=samps
  out
}


#####################################################
########## sampling of pi #########################
#######################################################

samp_pi=function(samp,Params,temper)
{
  L=samp$L
  temp=apply(L,1,function(x){all(x==2)})
  n=sum(temp)
  m=Params$J-n
  rbeta(1,(n+Params$a_pi+temper-1)/temper,(m+Params$b_pi+temper-1)/temper)
}


##################################################################
########## tune theta parameter ###############################
##################################################################
adaptive_tune=function(S,low=0.4,high=0.6)
{
  d1=dim(S)[1]
  d2=dim(S)[2]
  rej_rate=matrix(0,d1,d2)

  if( d2 > 1){ # for theta
    for(d in 1:d2)
    {
      temp=t(S[,1,])
      mcmcS=mcmc(temp)
      rej_rate[,d]=rejectionRate(mcmcS)
    }
    rej_rate=mean(rej_rate)
  } else { # for omega
    mcmcS = mcmc(S)
    rej_rate = rejectionRate(mcmcS)
  }

  if(rej_rate>high)
  {out=1} else if (rej_rate<low)
  {out=-1} else {out=0}
  out2=list(move=out,rate=rej_rate)
  out2
}


################################################################
##### calculate log-likelihood of a set of parameters ###########
####################################################################
log_dbeta_binom = function(X,D,P,over_dispersion ) 
{
  log_choose = lchoose(D,X) 
  numer = lbeta(X+over_dispersion*P,(D-X)+over_dispersion*(1-P))
  denom = lbeta(over_dispersion*P , over_dispersion*(1-P))
  
  log_choose + numer - denom 
}


prob_XD=function(X,D,Frac,L=NULL,Z,phi=NULL)
{
  M=L%*%Frac  # averaged copy number
  M2=Z%*%Frac  # averaged variant copy
  
  #### prob_X
  p=M2/M    # predicted VAF
  prob_X = dbinom(X,D,p,log = T)
  
  #### prob_D
  lamb=M%*%diag(as.vector(phi))/2
  prob_D=dpois(D,lamb,log=T)
  
  #### total probability 
  s=sum(prob_X)+sum(prob_D)
  
  s
}


prob_XD2=function(X,D,Frac,L=NULL,Z,phi=NULL) 
{
  M=L%*%Frac  # averaged copy number
  M2=Z%*%Frac  # averaged variant copy
  
  #### prob_X
  p=M2/M    # predicted VAF
  prob_X = dbinom(X,D,p,log = T)
  
  
  #### prob_D
  lamb=M%*%diag(as.vector(phi))/2
  prob_D=dpois(D,lamb,log=T)
  
  #### total probability 
  s=prob_X + prob_D
  
  s
}



