
################################################################
############# load packages ##################################
################################################################

require(tidyr)
require(copynumber)
require(ggplot2)
require(reshape)
require(dplyr)
require(coda)
require(gtools)
require(Rcpp)
require(RcppArmadillo)
require("igraph")


########################################################
######## load functions ################################
########################################################

source('./assist_fun.R')
source('./main_fun_II.R')
source('./tree_samp_fun.R')
source('./par_samp.R')
source('./call_seg.R')
sourceCpp("params.cpp")
source("Visualization.R")



#################################################################
########## MODEL INPUT ############################################
#################################################################

load("example.Rdata")
myseed = 1               # set random seed
foldername= "temp_out"         # set output foldername
dir.create(foldername)  # folder where outputs are saved

D=as.matrix(obs_data$D)  # total reads, J * T matrix
X=as.matrix(obs_data$X)  # variant reads, J * T matrix
mut_loc = obs_data$loc        # loci location matrix
segments = obs_data$segments  # segments



##############################################
######## load parameter file ################
############################################
phi=apply(D,2,median) # use median reads as read depth
phi=matrix(phi,nrow=1)

source("specify_pars.R")
par_disp(Params, MCMC_par)

##############################################
######## sampling ##########################
############################################
source("sampler.R")


##################################################
########## model selection:   ###############
##################################################

source("Model_select.R")
pdf(paste(foldername,"/","selection.pdf",sep=""))
BFE_calc(X,D,Temperature,foldername)
dev.off()


##################################################
########## visualization   ###############
##################################################
cat("Visualizing sampling results: \n")
Fit_visual(foldername,X,D)








