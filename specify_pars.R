################################################################################
#
# this file is used to specify model parameters
#
########################################################################################


########################################################
############## MODEL PARAMETERS ############################
########################################################

set.seed(myseed)

MCMC_par=Params=list()

#### number of samples
Params$T0= ncol(D)

#### number of genes
Params$J=nrow(D)

#### maximum number of copy
Params$max_CN=4

#### maximum number of mutant copies
Params$max_mut=2

##### hyper parameter of F: theta
Params$r=1.5

##### hyper parameter of pi
Params$a_pi=10^4
Params$b_pi=1

##### hyper parameter in prior of Z: prob of gaining one mutant allele
Params$zeta=10^-2

##### read depth parameter
Params$phi=phi

#### over_dispersion hyper-parameter
Params$lambda = 10


########################################################
############## get segments ############################
########################################################

if(is.null(segments)){
  cat("Genome segment information not found. Generating segments using loci position and total reads\n")
  if(is.null(mut_loc)){
    stop("Loci position information not found, cannot call segments")
  }
  segments=CallSeg(D,mut_loc,foldername)
}

Params$segments=segments


########################################################
############  specify MCMC parameters ##################
########################################################




 # MCMC_par$burnin=4000  # burnin sample size
 # MCMC_par$Nsamp=4000   # number of samples for inference
 # MCMC_par$Ntune=2000  # number of samples used for adaptive parameter tuning
 MCMC_par$burnin=150  # burnin sample size
 MCMC_par$Nsamp=150   # number of samples for inference
 MCMC_par$Ntune=150  # number of samples used for adaptive parameter tuning



MCMC_par$swap_interval=30 # make Matroplis Hastings move in every how many samples
MCMC_par$Nchain=5 # number of paralel chains
MCMC_par$delta_T=0.35

 Temperature=seq(1,by=MCMC_par$delta_T,length.out = MCMC_par$Nchain )  # temperatures

######################################################
############## adaptive tuning parameter ##############
######################################################

adapt_par=list(theta_tune=6)
adapt_par=lapply(c(1:MCMC_par$Nchain),function(x){adapt_par})



######################################################
##########  sequence of candidate Ks ##########
######################################################

# candidate subclone numbers K
# need K >= 2
Nclone=c(3:7) 
