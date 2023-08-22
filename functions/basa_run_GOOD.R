# basa_run.r
# Created by John Trochta
# Date created:  06/08/2019
# Summary:
# This script runs the Bayesian ASA model using the No-U-turn (NUTS) MCMC sampler to obtain posteriors.
# The adnuts package developed by Cole Monnahan to run NUTS with ADMB (hence adnuts) is used.
# I highly recommend users should read the following for more details on NUTS application to stock assessments:
#   Cole C Monnahan, Trevor A Branch, James T Thorson, Ian J Stewart, Cody S Szuwalski, 
#         Overcoming long Bayesian run times in integrated fisheries stock assessments, 
#         ICES Journal of Marine Science, Volume 76, Issue 6, November-December 2019, 
#         Pages 1477–1488, https://doi-org.eres.qnl.qa/10.1093/icesjms/fsz059

# Some annotated guidance on diagnostic checks is provided below:
# 1) If divergent transitions>0, increase target acceptance rate to reduce step size
#    e.g. control = list(adapt_delta = 0.9)
# 2) IF extreme global correlations, pass dense matrix estimated from previous run
#    e.g. control-list(metric=M) where M is matrix in untransformed space
#    for ADMB models, use MLE covairance with control=list(metric="mle")

# With this adnuts, most important diagnostics are the:
#   1) ESS (accounts for autocorrelation)-500 ESS is sufficient for most quantities
#   2) Potential Scale reduction (R hat)-R hat fails if >1.1
#   3) No max tree depths exceeded (<12)
#   4) <1% divergences
# reps <- parallel::detectCores()-1 # chains to run in parallel

#################################################################
library(data.table)

#################################################################
# Are you running this on a PC or a Mac
OS <- "MAC"

# BE SURE TO CHECK YOUR DIRECTORY
#template_files <- here::here("results/real_fits/base_ASA_2020_run_test/")  
template_files <- here::here("model/")
function_dir <- here::here("functions/")
setwd(template_files)

source(file=paste0(function_dir,"data_reader.R"))
source(file=paste0(function_dir,"data_header_reader.R"))

system("admb -s PWS_ASA")
# shell('admb PWS_ASA')


read.ADMB.files <- function(){
  # Parameters that are fixed and thus included as data
  #Z_3_8 <- 0.25
  egg_add <- 0.4
  
  # Store the file names from which data is available
  filename <- vector(length=2)
  filename[1]="PWS_ASA.dat"
  filename[2]="PWS_ASA(ESS).ctl"
  
  function_dir <- here::here("functions/") 
  source(file=paste0(function_dir,"data_reader.R"))
  source(file=paste0(function_dir,"data_header_reader.R"))
  
  PWS_ASA.dat <- data_reader(filename=filename[1]) # This is nyr - we want to start at nyr_tobefit
  PWS_ASA.dat <- PWS_ASA.dat[-1]
  nyr <- PWS_ASA.dat[[1]]
  nage <- PWS_ASA.dat[[2]]
  w_a_a <- PWS_ASA.dat[[3]] # From PWS_ASA.dat
  fecun <- PWS_ASA.dat[[4]]
  pc <- PWS_ASA.dat[[5]]
  pk <- PWS_ASA.dat[[6]]
  fbc <- PWS_ASA.dat[[7]]
  gc <- PWS_ASA.dat[[8]]
  sc <- PWS_ASA.dat[[9]]
  f_sp <- PWS_ASA.dat[[10]]
  
  mdm <- PWS_ASA.dat[[11]]
  egg <- PWS_ASA.dat[[12]]
  cv_egg <- PWS_ASA.dat[[13]]
  hydADFG <- PWS_ASA.dat[[15]]
  hydPWSSC <- PWS_ASA.dat[[17]]
  cv_hydPWSSC <- PWS_ASA.dat[[18]]
  seine <- PWS_ASA.dat[[19]]
  spac <- PWS_ASA.dat[[20]]
  juv <- PWS_ASA.dat[[21]]
  
  sero_comp <- PWS_ASA.dat[[23]]
  
  PWS_ASA_ESS.ctl <- data_reader(filename=filename[2])
  ESS_Se <- PWS_ASA_ESS.ctl[[1]]
  ESS_Sp <- PWS_ASA_ESS.ctl[[2]]
  ESS_Sero <- PWS_ASA_ESS.ctl[[3]]
  
  seine_indices <- which(rowSums(seine[1:nyr,])>0)
  spawnsurvey_indices <- which(rowSums(spac[1:nyr,])>0)
  egg_indices <- which(egg[1:nyr]>0)
  hydADFG_indices <- which(hydADFG[1:nyr]>0)
  hydPWSSC_indices <- which(hydPWSSC[1:nyr]>0)
  juvenile_indices <- which(juv[1:nyr]>0)
  
  sero_indices <- which(apply(sero_comp,2,function(x) all(x>=0)))
  
  model.data <- list(nyr=nyr,
                     nage=nage,
                     w_a_a=w_a_a[1:nyr,],
                     fecun=fecun[1:nyr,],
                     pc=pc[1:nyr,],
                     pk=pk,
                     fbc=fbc[1:nyr,],
                     gc=gc[1:nyr,],
                     sc=sc[1:nyr],
                     f_sp=f_sp[1:nyr],
                     ESS_Se=ESS_Se,
                     ESS_Sp=ESS_Sp,
                     ESS_Sero=ESS_Sero,
                     seine=seine[1:nyr,],
                     spac=spac[1:nyr,],
                     mdm=mdm[1:nyr],
                     egg=egg[1:nyr],
                     cv_egg=cv_egg[1:nyr],
                     hydADFG=hydADFG[1:nyr],
                     hydPWSSC=hydPWSSC[1:nyr],
                     cv_hydPWSSC=cv_hydPWSSC[1:nyr],
                     seine_indices=seine_indices,
                     spawnsurvey_indices=spawnsurvey_indices,
                     egg_indices=egg_indices,
                     hydADFG_indices=hydADFG_indices,
                     hydPWSSC_indices=hydPWSSC_indices,
                     #Z_3_8=Z_3_8,
                     egg_add=egg_add,
                     juv=juv[1:nyr],
                     juvenile_indices=juvenile_indices,
                     sero_comp=sero_comp[1:nyr,],
                     sero_indices=sero_indices)
  return(model.data)
}
model.data <- read.ADMB.files()
# Read in measured age comps
# The following commands skips over lines of the file to read in tables 
# so change the number skipped if file is modified  
nyr <- model.data$nyr # 35 for data up through 2014
seine.ac <- model.data$seine
spawn.ac <- model.data$spac
sero.ac <- model.data$sero_comp

# Read in the actual sample sizes
seine.SS <- read.table("agecomp_samp_sizes.txt",header=FALSE,skip=4,nrows=nyr)
spawn.SS <- read.table("agecomp_samp_sizes.txt",header=FALSE,skip=4+nyr+1,nrows=nyr)
sero.SS <- read.table("agecomp_samp_sizes.txt",header=FALSE,skip=4+nyr+1+nyr+1,nrows=nyr)

# Create empty matrices to fill estimated ESS and age comps
nage <- 10 # number of age classes (age 3 to 9+)

#Seine
Seine.ess.its <- matrix(0, nyr, 1) # Matrix to hold all iterations of the routine
Seine.ess.its <- seine.SS # fill in the first column with the recorded sample size 

#Spawn
Spawn.ess.its <- matrix(0, nyr, 1)
Spawn.ess.its <- spawn.SS

#Seroprevalence
Sero.ess.its <- matrix(0, nyr, 1)
Sero.ess.its <- sero.SS

# Change phases of the ESS in phases file to use the PWS_ASA(ESS_estimate)
ph <- readLines("PWS_ASA(phases).ctl",-1)
ph[4] <- 1
writeLines(ph,"PWS_ASA(phases).ctl")

# LOOP THROUGH AND ITERATIVELY CALCULATE ESS
convergence <- 0

seine.ESS <- seine.SS
spawn.ESS <- spawn.SS
sero.ESS <- sero.SS
its <- 1

#while(convergence==0) {
for(i in 1:2){
  # Create "PWS_ASA(ESS_estimate).ctl" with sample sizes (the original sample sizes on the first iteration)
  write.table(rbind("# PWS age comp effective sample sizes","# Seine ESS", seine.ESS,
                    " ", "# Spawn ESS", spawn.ESS,
                    " ", "# Sero ESS", sero.ESS),
              file = "PWS_ASA(ESS_estimate).ctl", append = F, sep = " ",
              row.names=FALSE,col.names=FALSE,quote=F)
  
  # Compile and Run PWS_ASA
  if(OS=="MAC"){
    system("./PWS_ASA -pinwrite -nohess")
  }else if(OS=="PC"){
    shell('PWS_ASA  -pinwrite -nohess')
  }
  
  
  # Read in the estimated seine and spawner age comps
  seine.ac.est <- read.table("SeAC_pd.rep", header = FALSE) 
  spawn.ac.est <- read.table("SpAC_pd.rep", header = FALSE)
  sero.ac.est <- read.table("SeroAC_pd.rep", header = FALSE)
  
  # Calculate the ESS
  seine.ESS <- rowSums(seine.ac.est*(1-seine.ac.est))/rowSums((seine.ac-seine.ac.est)^2)
  spawn.ESS <- rowSums(spawn.ac.est*(1-spawn.ac.est))/rowSums((spawn.ac-spawn.ac.est)^2)
  sero.ESS <- rowSums(sero.ac.est*(1-sero.ac.est))/rowSums((sero.ac-sero.ac.est)^2)
  
  # Remove the missing years of age comps
  seine.ESS.rem <- seine.ESS[!(seine.ac[,1]==-9)]
  spawn.ESS.rem <- spawn.ESS[!(spawn.ac[,1]==-9)]
  sero.ESS.rem <- sero.ESS[!(sero.ac[,1]==-9)]
  
  # Calculate the ratio of ESS to original sample sizes
  seine.ratio <- seine.ESS.rem/seine.SS[seine.ac[,1]!=-9,1]
  spawn.ratio <- spawn.ESS.rem/spawn.SS[spawn.ac[,1]!=-9,1]
  sero.ratio <- sero.ESS.rem/sero.SS[sero.ac[,1]!=-9,1]
  
  # Calculate the harmonic means
  seine.hm <- 1/mean(1/seine.ratio)
  spawn.hm <- 1/mean(1/spawn.ratio)
  sero.hm <- 1/mean(1/sero.ratio)
  
  # Compare this harmonic mean to the previous using a convergence criteria (WHAT AM I CONVERGING!!!!)
  if(its==1) {
    convergence <- 0
    seine.hmS <- seine.hm
    spawn.hmS <- spawn.hm
    sero.hmS <- sero.hm
  }else{
    seine.test <- abs(seine.hm - seine.hmS[its-1])/seine.hmS[its-1]*100
    spawn.test <- abs(spawn.hm - spawn.hmS[its-1])/spawn.hmS[its-1]*100
    sero.test <- abs(sero.hm - sero.hmS[its-1])/sero.hmS[its-1]*100
    convergence <- (seine.test<.1 & spawn.test<.1 & sero.test<0.1) # This criteria was arbitrarily chosen (0.1% change)
    seine.hmS <- rbind(seine.hmS,seine.hm)
    spawn.hmS <- rbind(spawn.hmS,spawn.hm) 
    sero.hmS <- rbind(sero.hmS,sero.hm) 
  }
  
  # Now multiply the harmonic mean by the sample size to get the new ESS 
  seine.ESS <- round(seine.hm*seine.SS, digits=0)
  spawn.ESS <- round(spawn.hm*spawn.SS, 0)
  sero.ESS <- round(sero.hm*sero.SS, 0)
  
  # Use the average ESS for all years (each years obs weighted equally)
  # seine.ESS[seine.ESS>0] <- round(mean(seine.ESS[seine.ESS>0]), digits=0)
  # spawn.ESS[spawn.ESS>0] <- round(mean(spawn.ESS[spawn.ESS>0]), digits=0)
  
  # Denote the missing values
  seine.ESS[(seine.ac[,1]==-9),1] <- -9
  spawn.ESS[(spawn.ac[,1]==-9),1] <- -9
  sero.ESS[(sero.ac[,1]==-9),1] <- -9
  #     seine.ESS <- seine.SS
  #     spawn.ESS <- spawn.SS
  #     seine.ESS[(seine.ac[,1]==-9)] <- 0
  #     spawn.ESS[(spawn.ac[,1]==-9)] <- 0
  
  sero.ESS <- sero.SS
  
  # Fill in this iteration's ESS
  Seine.ess.its <- cbind(Seine.ess.its,round(seine.ESS,0))
  Spawn.ess.its <- cbind(Spawn.ess.its,round(spawn.ESS,0))
  Sero.ess.its <- cbind(Sero.ess.its,round(sero.ESS,0))
  
  # Cease iterations if convergence hasn't happened after so many...
  if(its==10) {
    break
  }
  its <- its+1
}

# Turn of the phases for the ESS calculation so it no longer recalculates ESS in future model runs
ph <- readLines("PWS_ASA(phases).ctl",-1)
ph[4] <- -1
writeLines(ph,"PWS_ASA(phases).ctl")

# Now write the converged ESS to a ctl file to be used for model runs
write.table(rbind("# PWS age comp effective sample sizes",paste0("# (",date(),")")," ",
                  "# Seine ESS", seine.ESS,
                  " ", "# Spawn ESS", spawn.ESS,
                  " ", "# Sero ESS", sero.ESS),
            file = "PWS_ASA(ESS).ctl", append = F, sep = " ",
            row.names=FALSE,col.names=FALSE,quote=F)
  
######################################################
library(adnuts)
library(snowfall)
library(rstan)
library(r4ss)

setwd(template_files)
reps <- 4
set.seed(8558)
seeds <- sample(1:1e4, size=reps)
system('./PWS_ASA -pinwrite')
#system('./PWS_ASA -pinwrite -hbf 1 -nox -iprint 200 -mcmc 15')

######################################################
# Random Walk metropolis for very initial exploration

# Need the Hessian
# if(OS=="MAC"){
#   system("admb -s PWS_ASA")
#   system("./PWS_ASA")
# }else if(OS=="PC"){
#   shell('admb PWS_ASA')
#   shell('PWS_ASA -hbf 1')
# }

# thin <- 5000
# iter <- 1000*thin
# warmup <- iter/5
# inits <- NULL ## start chains from MLE
# pilot <- sample_admb("./PWS_ASA", iter=iter, thin=thin, seeds=seeds, init=inits,
#                      parallel=TRUE, chains=reps, warmup=warmup,
#                      cores=reps, algorithm='RWM')
# launch_shinyadmb(pilot)
# ## Check convergence
# mon <- monitor(pilot$samples, warmup=pilot$warmup, print=FALSE)
# max(mon[,'Rhat'])
# min(mon[,'n_eff'])
# ## Examine the slowest mixing parameters
# slow <- names(sort(mon[,'n_eff']))[1:8]
# pairs_admb(fit=pilot, pars=slow)
######################################################
# NUTS starting from MLE
# setwd('..')
# 
# nuts.mle <-
#   sample_admb(model="./PWS_ASA", iter=500, init=NULL, algorithm='NUTS',  seeds=seeds,
#               parallel=TRUE, chains=reps, warmup=100, cores=reps,
#               control=list(metric="mle", adapt_delta=0.99))
# launch_shinyadmb(nuts.mle)
# saveRDS(fit, file=paste0("results/NUTS_", m, ".RDS"))
######################################################
pars <- read.admbFit('PWS_ASA')
inits <- list()
for(j in 1:reps){
  # This is a check on the CV of pars - basically high uncertainty indicates parameter is fixed or uninformed
  # attempts to avoid boundary issues in further estimation
  if(pars$std[1]/pars$est[1]<100 | !is.infinite(pars$std[1]/pars$est[1])){
    inits[[j]] <- rnorm(1,pars$est[1],sd=pars$std[1])  
  }else{inits[[j]] <- pars$est[1]}
  
  for(k in 2:length(pars$est)){
    if(pars$names[k]=="annual_inf"){
      inits[[j]] <- c(inits[[j]],runif(1,0,1))  
    }else if(pars$names[k]=="recov_prob"){
      inits[[j]] <- c(inits[[j]],runif(1,0.2,0.9))  
    }else if(pars$std[k]/pars$est[k]<100 | !is.infinite(pars$std[k]/pars$est[k])){
      inits[[j]] <- c(inits[[j]],rnorm(1,pars$est[k],sd=pars$std[k]))  
    }else{inits[[j]] <- c(inits[[j]],pars$est[k])}
  }
  
  inits[[j]] <- inits[[j]][-length(inits[[j]])]
}

# ADMB command for running nuts
# PWS_ASA -nox -noest -nohess -maxfn 0 -nuts -mcmc 2000 -warmup 500 -chain 1 -mcseed 8682524 -max_treedepth 12 -adapt_delta 0.8 -adapt_mass -mcpin init.pin

# Pilot run to check 
start_time <- Sys.time()
fit.1 <- sample_nuts(model='./PWS_ASA',path=template_files,
                       iter=4000, init=inits, seeds=seeds,
                       chains=reps,parallel=TRUE,cores=reps,
                       warmup=500, mceval=TRUE, control=list(adapt_delta=0.925))
# pilot <- sample_admb(model='./PWS_ASA',path=template_files,
#                    iter=500, init=inits, algorithm='NUTS',  seeds=seeds,
#                    parallel=TRUE, chains=reps, warmup=100, cores=reps,
#                    mceval=TRUE, control=list(adapt_delta=0.9))
end_time <- Sys.time()
end_time - start_time

# 11/17/2020:  36.58 minutes
# 11/17/2020:  40.33 minutes for 2020 data
# 12/04/2020:  42.61 minutes for basa_ASA_2020_run_2
# 12/04/2020:  40.59 minutes for basa_ASA_2020_run_3
# 12/08/2020:  43.58 minutes for basa_ASA_2020_run
# 12/08/2020:  44.70 minutes for base_ASA_2020_run_inf vul est
# 12/09/2020:  46.31 minutes for base_ASA_2020_run_ess50
# 12/09/2020:  45.03 minutes for base_ASA_2020_run_2_ess50
# 12/09/2020:  43.57 minutes for base_ASA_2020_run_3_ess50
# 01/01/2021:  43.09 minutes for base_ASA_2020_run_binomial_test

# Do my own monitoring
mon <- monitor(fit.1$samples, warmup=fit.1$warmup, print=FALSE)
x <- extract_sampler_params(fit.1)
data.frame(model=fit.1$model, pct.divs=100*sum(x$divergent__)/nrow(x))

write.csv(mon, file='table_parameter_posteriors.csv')

divergences <- sum(fit.1$sampler_params[[1]][,5])/length(fit.1$sampler_params[[1]][,5])
## Examine the slowest mixing parameters
# slow <- names(sort(mon[,'n_eff']))[1:8]
# pairs_admb(fit=fit, pars=slow)

sum(x$divergent__)/nrow(x)<=0.01
max(mon[,'Rhat'])<=1.1

sum_dia <- data.frame(divergences.from.extract.function=sum(x$divergent__)/nrow(x),
                      divergences.from.fit.list=sum(fit.1$sampler_params[[1]][,5])/length(fit.1$sampler_params[[1]][,5]),
                      min.ESS=min(mon[,'n_eff']),
                      which.min.ESS=names(which.min(mon[,'n_eff'])),
                      max.Rhat=max(mon[,'Rhat']),
                      which.max.Rhat=names(which.max(mon[,'Rhat'])),
                      time.elapsed=end_time-start_time)
write.table(sum_dia,file='table_convergence_diagnostics.csv',sep=",",row.names=FALSE)
saveRDS(fit.1, file="NUTS_fit.RDS")

# Launch Shiny App to check diagnostics online
launch_shinyadmb(fit.1)

rm(list = ls(all.names = TRUE))

