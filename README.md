_SIFA_ (tumor **S**ubclone **I**dentification by **F**eature **A**llocation) is a Baysian method to identify tumor subclones using WGS data. This page will guide you through the basic steps of using _SIFA_.

Currently SIFA requires sample size to be at least two, since a unique tree cannot be identified with only one sample.

## Software dependencies
_SIFA_ is written in `R` and `C++`. Please install the following packages in R prior to implementing our software:

- data manipulation: `tidyr`,`reshape`,`dplyr`
- segment calling: `copynumber`
- Bayesian analysis: `coda`
- integrating `c++` functionality: `Rcpp`,`RcppArmadillo`
- visualization: `ggplot2`, `igraph`
- others: `gtools`


## Prepare software input

_SIFA_  takes an `.Rdata` file as input. 

The `.Rdata` file should contain a list `obs_data` with required fields: `obs_data$D` for total reads matrix, `obs_data$X` for mutant reads matrix, and `obs_data$loc` for mutation location matrix, and one optional field: `obs_data$segments` for loci segmentation matrix. 

The total reads (`D`) and mutation location (`loc`) data will be used to calculate genome segments, unless the optional input `obs_data$segments` is given.

**Required inputs**:

- `obs_data$D` and `obs_data$X` should take the following format:

	loci | sample 1 | sample 2 | sample 3 | ...
  ------- | --------- | -------- | -------- | --------
  locus 1 | 12 | 11 | 33  | ...
   locus 2 | 5 | 8 | 7 | ...
  ...    | ... | ... | ... | ... 
  locus J | 22 | 10 | 17 | ... 

- `obs_data$loc` should take the following format:

   chromosome | position | gene 
  ------- | --------- | -------- 
  1 | 7660469 |  CAMTA1
  3 | 88482840 | 
  13 | 102703724  |  FGF14 
  ...    | ... | ... 
  23 | 153383479 | 
For loci in non-coding regions, the `gene` column can be left blank. 

**Optional inputs:**

- `obs_data$segments` should take the following format:

  segments | start | end 
  ------- | --------- | -------- 
  segment 1 | 1 |  5 
  segment 2 | 6 |  25
  segment 3 | 26  |  40
  ...    | ... | ... 
  segment S |  155 | J
Each row of the matrix represents one segment, with the two entries marking the starting and ending locus of the segment.

## Using _SIFA_

To use _SIFA_, please set R working directory to `SIFA` after cloning this repository. Make sure you have all the dependencies correctly installed, then run source code `SIFA_app.R` line by line.

- In the _MODEL INPUT_ section of the code, load the `.Rdata` where your inputs are located, specify random seed `myseed`, and specify the folder `foldername` to store output files (a new folder will be created if it does not exist). For example: 

	```r
#############################################
########## MODEL INPUT ######################
#############################################
load("example.Rdata")
myseed = 1                # set random seed
foldername = "temp_out"   # set output foldername
	```

- Next, you need to specify Bayesian sampling parameters in `specify_pars.R`. For most of the parameters, default values work just fine. Some of the parameters you can change are:

	```r
	#### maximum number of copy
	Params$max_CN=4
	#### maximum number of mutant copies
	Params$max_mut=2
	
	#### MCMC sampling parameters 
	MCMC_par$burnin=4000  # burnin sample size
	MCMC_par$Nsamp=4000   # number of samples for inference
	MCMC_par$Ntune=2000  # number of samples used for 	adaptive parameter tuning
	Nclone=c(3:7) # candidate subclone numbers K
	```

- Implement the following files:
	- `sampler.R` to perform sampling
	- `DIC_select.R` to perform model selection
	- `Visualization.R` for results visualization   
	Visualization results will list top 3 frequent trees (when more than 3 tree structures exist) in posterior samples, and display corresponding parameter estimations.

During the sampling process, all samples for each individual K will be stored in one `.Rdata` file.
## Contact
Please feel free to contact <li.zeng@yale.edu> if you have any question.