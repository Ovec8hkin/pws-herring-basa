# Script for performing retrospective analyses on BASA
# Joshua Zahner | 05/17/2022 
#
# Creates a new subdirectory ('retrospectives/') and runs
# BASA on successivley fewer years of data (up to a provided
# number of 'peels'). Final biomass estimates and recrtuiment
# deviates are aggregated and analysed for bias patterns relative
# to the most recent (current year) full BASA run.
library(here)
library(doParallel)
library(tidyverse)
library(icesAdvice)
setwd(here::here())
source(file=file.path(here::here(), "functions/fun_read_dat.R"))
source(file=file.path(here::here(), "functions/fun_write_dat.R"))
source(file=file.path(here::here(), "functions/run_basa.R"))

n.peels <- 5
main.basa.directory <- paste0(here::here("model/"))

run.parallel <- FALSE
total.cores <- parallel::detectCores()
parallel.runs <- (total.cores-1) %/% 4      # BASA runs using 4 nodes

run.retropspective <- function(i, ...){
    new.year <- 2022-i
    new.dir.name <- paste0("basa_", new.year)
    if(!dir.exists(new.dir.name)){
        dir.create(new.dir.name)
    }
    setwd(new.dir.name)
    
    if(!file.exists("mcmc_out/PFRBiomass.csv")){
        print(paste("mcmc_out directory does not exist in", new.dir.name))
        file.symlink(paste0(main.basa.directory, "PWS_ASA.TPL"), ".")
        file.symlink(paste0(main.basa.directory, "PWS_ASA(par).ctl"), ".")
        file.symlink(paste0(main.basa.directory, "PWS_ASA(phases).ctl"), ".")
        file.symlink(paste0(main.basa.directory, "PWS_ASA(sim_settings).ctl"), ".")
        dat.files <- read.data.files(main.basa.directory)
        fun_write_dat(dat.files, i)
        run.basa(paste0("retrospectives/", new.dir.name))    
    }else{
        print("mcmc_out directory already exists")
    }
    setwd("..")
}

create.retrospective.matrix <- function(n.peels, nyr, cyr=2023, fname="PFRBiomass.csv"){
    # Aggregate biomass estimates for each retrospective run
    annual.estimate <- matrix(NA, nyr, n.peels+1)
    for(i in n.peels:1){
        print(i)
        est <- read_csv(paste0(here::here(), "/retrospectives/", "basa_", cyr-i,"/mcmc_out/", fname), col_names=FALSE, show_col_types = FALSE) %>%
                summarise(
                    across(
                    everything(),
                    median
                    )
                ) 
        
        est <- t(as.matrix(est))
        rownames(est) <- NULL
        for(j in 1:length(est)){
            print(j)
            annual.estimate[j, i+1] <- est[j, 1]
        }
    }

    # Special case to get the current assessment run which lives in a different directory
    annual.estimate[1:nyr,1] <- as.numeric(read_csv(paste0(here::here(), "/model/mcmc_out/", fname), col_names=FALSE, show_col_types = FALSE) %>%
                                            summarise(
                                                across(
                                                everything(),
                                                median
                                                )
                                            ))
    
    return(annual.estimate)
}

# Set up retrospectives directory
if(!dir.exists("retrospectives/")){
    dir.create("retrospectives/")
}
setwd("retrospectives/")

# Need to copy and then modify all of input files so they only reflect the data
# that was available in that year. Then actually run the assessment on all of
# the old datasets. Do this in parallel to speed things up.
if(run.parallel){
    cluster <- parallel::makeCluster(parallel.runs, type="FORK")
    registerDoParallel(cluster)
    foreach(i=1:n.peels) %dopar% {
        run.retropspective(i)
    }
    stopCluster(cluster)
}else{
    for(i in 1:5){
        run.retropspective(i=i)
    }
}

########################## Analaysis and Plotting Code ###########################

annual.biomass.estimate <- create.retrospective.matrix(4, 45, cyr=2022)

mohns.rho <- icesAdvice::mohn(annual.biomass.estimate, details=TRUE)
rel.bias <- c(mohns.rho$compare[,"relbias"], NA)
rho <- mohns.rho$rho

# Compute the peel estimates for each run
final.year.est <- annual.biomass.estimate[38:nrow(annual.biomass.estimate), 1]
peel.estimates <- apply(as.matrix(annual.biomass.estimate), 2, function(x) x[max(which(!is.na(x)))])
peel.estimates <- rev(peel.estimates)


# Plot retrospective plot (can also be donee with icesAdvice::mohn(..., plot=TRUE))
rownames(annual.biomass.estimate) <- 1980:2024
colors <- palette(rainbow(n.peels+1))

plot(1980:2024, rep(NA, 45), ylim=c(0, 150000), xlab="Year", ylab="SSB (mt)", main="Retrospective Analaysis of SSB")
for(p in 1:n.peels){
    data <- annual.biomass.estimate[,1+p]
    lines(1980:2024, data, col=colors[p], lwd=2)
    points(2024-p, rev(peel.estimates)[1+p], col=colors[p], pch=19, lwd=10)
}
lines(1980:2024, annual.biomass.estimate[,1], col="black", lwd=2.3)
text(2019.5, 110000, paste0("Mohn's Rho: ", round(rho, 3)))
abline(h=19958, lty=2)
abline(h=39885, lty=2)

# Fancy thing for adding the relative biases to the legend
legs <- rep(NA, length(2018:2022))
for(i in 1:length(legs)){
    legs[i] <- paste0(2022-n.peels-1+i, " (bias = ", round(rel.bias[i], 3), ")")
}
legend(x=2014.5, y=150000, legend=rev(legs), col=c("black", colors), lwd=2)

# Make squid plot
plot(1:n.peels, rep(NA, n.peels), xlim=c(0, n.peels+1), ylim=c(0, 150000), xlab="")
colors <- palette(rainbow(n.peels+1))
for(i in 1:(n.peels+1)){
    y <- 2022-1-n.peels+i
    data <- rev(annual.biomass.estimate[as.character(y),])
    non.na.index <- which(!is.na(data))[1]
    data <- data[non.na.index:length(data)]
    lines(0:length(data), c(0, data), col=colors[i], lwd=3)
    text(x=length(data), y=data[length(data)]+2500, y)
}

# Make reverse squid plot
plot(1:n.peels, rep(NA, n.peels), xlim=c(0.5, n.peels+1), ylim=c(-2, 2), xlab="")
colors <- palette(rainbow(n.peels+1))
for(i in 1:(n.peels)){
    y <- 2022-1-n.peels+i
    data <- rev(annual.biomass.estimate[as.character(y),])
    data <- (data - data[length(data)])/data[length(data)]
    non.na.index <- which(!is.na(data))[1]
    data <- data[non.na.index:length(data)]
    lines(1:length(data), c(data), col=colors[i], lwd=3)
    text(x=0.75, y=data[1], y)
}
abline(h=0, lty=2)
