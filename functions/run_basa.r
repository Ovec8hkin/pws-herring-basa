# Sets up estimated sample sizes and runs BASA using the NUTS
# sample in ADMB (through adnuts package).
#
# Joshua Zahner | 05/17/2022
library(data.table)
library(tidyverse)
library(adnuts)
library(snowfall)
library(rstan)
library(r4ss)

function.dir <- here::here("functions/")
source(file=paste0(function.dir, "/data_reader.R"))
source(file=paste0(function.dir, "/data_header_reader.R"))
source(file=paste0(function.dir, "/read_admb_dat_files.R"))
source(file=paste0(function.dir, "/calculate_ess.R"))
source(file=paste0(function.dir, "/init_admb_params.R"))

#################################################################
# Are you running this on a PC or a Mac
OS <- "MAC"

run.basa <- function(model.dir, n.samples=2000, n.warmup=700, n.time=5){
    # BE SURE TO CHECK YOUR DIRECTORY
    template.files <- here::here(model.dir)
    print(template.files)
    setwd(template.files)
    system("admb -s PWS_ASA")
    #shell("admb -s PWS_ASA")

    model.data <- read.admb.files()

    nyr.fit <- model.data$nyr.fit # 44 for data up through 2023
    nyr.tot <- model.data$nyr.tot # 44 for data up through 2023
    nage <- 10                    # number of age classes (age 0 to 9+) (only care about 3-9+)

    # Read in measured age comps
    seine.age.comp <- model.data$seine.age.comp
    spawn.age.comp <- model.data$spawn.age.comp
    sero.age.comp <- model.data$sero.comp

    # Read in the actual sample sizes
    # The following commands skips over lines of the file to read in tables 
    # so change the number skipped if file is modified  
    seine.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4,                 nrows=nyr.tot)
    spawn.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+1*(nyr.tot+1),   nrows=nyr.tot)
    sero.samp.size  <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+2*(nyr.tot+1),   nrows=nyr.tot)


    # Create empty matrices to fill estimated ESS and age comps
    seine.ess.its <- matrix(0, nyr.tot, 1)  # Matrix to hold all iterations of the routine
    seine.ess.its <- seine.samp.size        # fill in the first column with the recorded sample size 

    spawn.ess.its <- matrix(0, nyr.tot, 1)
    spawn.ess.its <- spawn.samp.size

    sero.ess.its  <- matrix(0, nyr.tot, 1)
    sero.ess.its  <- sero.samp.size

    # Change phases of the ESS in phases file to use the PWS_ASA(ESS_estimate)
    phases <- readLines("PWS_ASA(ESS).ctl", -1)
    phases[5] <- 1
    writeLines(phases, "PWS_ASA(ESS).ctl")

    # LOOP THROUGH AND ITERATIVELY CALCULATE ESS
    convergence <- 0

    seine.ess <- seine.samp.size
    spawn.ess <- spawn.samp.size
    sero.ess  <- sero.samp.size
    its <- 1

    #age.comps <- list(seine=seine.age.comp, spawn=spawn.age.comp, vshv=vhsv.age.comp, ich=ich.age.comp)
    #start.ess <- list(seine=seine.ess, spawn=spawn.ess, vhvs=vhsv.ess, ich=ich.ess)
    #samp.size <- list(seine=seine.samp.size, spawn=spawn.samp.size, vhsv=vhsv.samp.size, ich=ich.samp.size)

    # calc.ess <- calculate.ess(age.comps, start.ess, samp.size, nyr.fit)

    for(i in 1:2){
      # Create "PWS_ASA(ESS_estimate).ctl" with sample sizes (the original sample sizes on the first iteration)
      write.table(rbind("# PWS age comp effective sample sizes","# Seine ESS", seine.ess,
                        " ", "# Spawn ESS", spawn.ess,
                        " ", "# Sero ESS", sero.ess),
                  file = "PWS_ASA(ESS_estimate).ctl", append = F, sep = " ",
                  row.names=FALSE,col.names=FALSE,quote=F)
      
      # Compile and Run PWS_ASA
      if(OS=="MAC"){
        system("./PWS_ASA -pinwrite -nohess")
      }else if(OS=="PC"){
        shell('PWS_ASA  -pinwrite -nohess')
      }
      
      
      # Read in the estimated seine and spawner age comps
      seine.age.comp.est <- read.table("SeAC_pd.rep", header = FALSE) 
      spawn.age.comp.est <- read.table("SpAC_pd.rep", header = FALSE)
      sero.age.comp.est <- read.table("SeroAC_pd.rep", header = FALSE)
      
      # Calculate the ESS
      seine.ess <- rowSums(seine.age.comp.est*(1-seine.age.comp.est))/rowSums((seine.age.comp-seine.age.comp.est)^2)
      spawn.ess <- rowSums(spawn.age.comp.est*(1-spawn.age.comp.est))/rowSums((spawn.age.comp-spawn.age.comp.est)^2)
      sero.ess <- rowSums(sero.age.comp.est*(1-sero.age.comp.est))/rowSums((sero.age.comp-sero.age.comp.est)^2)
      
      # Remove the missing years of age comps
      seine.ess.rem <- seine.ess[!(seine.age.comp[,1]==-9)]
      spawn.ess.rem <- spawn.ess[!(spawn.age.comp[,1]==-9)]
      sero.ess.rem <- sero.ess[!(sero.age.comp[,1]==-9)]
      
      # Calculate the ratio of ESS to original sample sizes
      seine.ratio <- seine.ess.rem/seine.samp.size[seine.age.comp[,1]!=-9,1]
      spawn.ratio <- spawn.ess.rem/spawn.samp.size[spawn.age.comp[,1]!=-9,1]
      sero.ratio <- sero.ess.rem/sero.samp.size[sero.age.comp[,1]!=-9,1]
      
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
      seine.ess <- round(seine.hm*seine.samp.size, digits=0)
      spawn.ess <- round(spawn.hm*spawn.samp.size, 0)
      sero.ess <- round(sero.hm*sero.samp.size, 0)
      
      # Use the average ESS for all years (each years obs weighted equally)
      # seine.ess[seine.ess>0] <- round(mean(seine.ess[seine.ess>0]), digits=0)
      # spawn.ess[spawn.ess>0] <- round(mean(spawn.ess[spawn.ess>0]), digits=0)
      
      # Denote the missing values
      seine.ess[(seine.age.comp[,1]==-9),1] <- -9
      spawn.ess[(spawn.age.comp[,1]==-9),1] <- -9
      sero.ess[(sero.age.comp[,1]==-9),1] <- -9
      #     seine.ess <- seine.samp.size
      #     spawn.ess <- spawn.samp.size
      #     seine.ess[(seine.age.comp[,1]==-9)] <- 0
      #     spawn.ess[(spawn.age.comp[,1]==-9)] <- 0
      
      sero.ess <- sero.samp.size
      
      # Fill in this iteration's ESS
      seine.ess.its <- cbind(seine.ess.its,round(seine.ess,0))
      spawn.ess.its <- cbind(spawn.ess.its,round(spawn.ess,0))
      sero.ess.its <- cbind(sero.ess.its,round(sero.ess,0))
      
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

    # Turn of the phases for the ESS calculation so it no longer recalculates ESS in future model runs
    # Now write the converged ESS to a ctl file to be used for model runs
    write.table(rbind("# PWS age comp effective sample sizes",paste0("# (",date(),")")," ",
                      "# Seine ESS", seine.ess,
                      " ", "# Spawn ESS", spawn.ess,
                      " ", "# Sero ESS", sero.ess),
                file = "PWS_ASA(ESS).ctl", append = F, sep = " ",
                row.names=FALSE,col.names=FALSE,quote=F)
      
    ######################################################
    # Create reps x starting par vectors, and run NUTS
    setwd(template.files)
    reps <- 4
    set.seed(8558)
    seeds <- sample(1:1e4, size=reps)
    #system("admb -s PWS_ASA")
    system("./PWS_ASA -pinwrite")

    inits <- init.admb.params(reps)

    # ADMB command for running nuts
    # PWS_ASA -nox -noest -nohess -maxfn 0 -nuts -mcmc 2000 -warmup 500 -chain 1 -mcseed 8682524 -max_treedepth 12 -adapt_delta 0.8 -adapt_mass -mcpin init.pin

    # Pilot run to check 

    if(dir.exists("mcmc_out")){
        system("rm mcmc_out/*.csv")
    }

    start.time <- Sys.time()
    fit.1 <- sample_nuts(model='./PWS_ASA',path=template.files,
                            iter=n.samples,
                            warmup=n.warmup,
                            #warmup=100,
                            duration = n.time,
                            init=inits, seeds=seeds, chains=reps,cores=reps,
                            mceval=TRUE,
                            control=list(
                                adapt_delta=0.9,
                                #max_treedepth=16,
                                metric="mle"
                            )
                        )
    end.time <- Sys.time()
    total.time <- end.time - start.time

    return(list(fit1=fit.1, time=total.time))
}