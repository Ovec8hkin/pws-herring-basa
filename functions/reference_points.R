source(file = file.path(here::here(), "functions", "fun_read_dat.R"))

dat.fname <- file.path(here::here(), "model", "PWS_ASA.dat")

M=0.20

dat <- data_reader(dat.fname)
waa <- dat[[4]]
waa <- apply(waa, 2, mean)
maa <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
vaa <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1)

pars <- read.par.file(file.path(here::here(), "model", "PWS_ASA.par"))

base.R0 <- pars$log_MeanAge0+mean(pars$annual_age0devs)

log.R0.high <- pars$log_MeanAge0+mean(pars$annual_age0devs[1:12])
log.R0.low <- pars$log_MeanAge0+mean(pars$annual_age0devs[13:32])
#R0 <- 1000000*exp(log.R0)

naa.per.recruit <- function(u=0.0, v){

    s <- exp(-M)

    n.max <- length(waa)
    naa <- rep(0, n.max)
    naa[1] <- 1
    for(a in 2:(n.max-1)){
        naa[a] <- naa[a-1]*(1-u*v[a])*s
    }
    naa[n.max] <- naa[a-1]*((1-u*v[a])*s)/(1-((1-u*v[a])*s))
    
    return(naa)
  
}

yield.per.recruit <- function(u=0.0, w, m, v){
  naa <- naa.per.recruit(u, v)
  YPR <- sum(naa*w*v*u)
  return(YPR)
}

spawn.biomass.per.recruit <- function(u=0.0, w, m, v){
  naa <- naa.per.recruit(u, v)
  f <- w*m
  SBPR <- sum(naa*f)
  return(SBPR)
}

spawning.potential.ratio <- function(u=0.0, w, m, v){
  sbpr <- spawn.biomass.per.recruit(u, w, m,v)
  sbpr0 <- spawn.biomass.per.recruit(0, w, m, v)
  spr <- sbpr/sbpr0
  return(spr)
}

npr <- naa.per.recruit(v=vaa)
ypr <- yield.per.recruit(w=waa, m=maa, v=vaa)
sbpr <- spawn.biomass.per.recruit(w=waa, m=maa, v=vaa)
spr <- spawning.potential.ratio(w=waa, m=maa, v=vaa)

unfished.spawning.biomass <- (sbpr*exp(log.R0.low))

thresholds <- c(0.40, 0.20)

unfished.spawning.biomass*thresholds

fishing.mortality <- seq(0.0, 0.50, 0.001)

yprs <- sapply(fishing.mortality, yield.per.recruit, w=waa, m=maa, v=vaa)
sbprs <- sapply(fishing.mortality, spawn.biomass.per.recruit, w=waa, m=maa, v=vaa)
sprs <- sapply(fishing.mortality, spawning.potential.ratio, w=waa, m=maa, v=vaa)

AUBs <- sapply(c(1:44), function(t) exp(pars$log_MeanAge0+mean(pars$annual_age0devs[1:t]))*spawn.biomass.per.recruit(w=waa, m=maa, v=vaa))



plot(x=fishing.mortality, yprs*exp(base.R0), type="l", col="red", ylab="Yield (mt)", xlab="F", main="Yield vs Fishing Mortality")
plot(x=fishing.mortality, sbprs*exp(base.R0), type="l", col="red", ylab="Spawning Biomass (mt)", xlab="F", main="Spawning Biomass vs Fishing Mortality")
plot(x=fishing.mortality, sprs, type="l", col="red", ylab="SPR", xlab="F", main="SPR vs Fishing Mortality")

plot(1980:2023, AUBs, type="l", ylab="Unfished Spawning Biomass (mt)", xlab="F", main="Unfished Spawning Biomass vs Fishing Mortality")

ypr2 <- function(u, w, m, v){
    return(-1*yield.per.recruit(u, w, m, v))
}

# 0.21336
optim(0.2, ypr2, w=waa, m=maa, v=vaa, lower=0.0, upper=1.0, method="Brent")

yield.per.recruit(0.21336, waa, maa, vaa)*exp(base.R0)
spawn.biomass.per.recruit(0.21336, waa, maa, vaa)*exp(base.R0)
spawning.potential.ratio(0.21336, waa, maa, vaa)

biomass.ests <- apply(read.biomass.estimates(here::here("model")), 2, median)

plot(1980:2023, AUBs, type="l", ylab="Unfished Spawning Biomass (mt)", ylim=c(0, 500000), xlab="year", main="Unfished Spawning Biomass vs Fishing Mortality")
lines(1980:2023, biomass.ests)

AUBs

data <- data.frame(year=1980:2023, AUB=AUBs, biomass=biomass.ests)
data$dep <- data$biomass/data$AUB
data$stable_AUB <- sapply(c(1:44), function(t) exp(pars$log_MeanAge0+mean(pars$annual_age0devs))*spawn.biomass.per.recruit(w=waa, m=maa, v=vaa))
data$stable_dep <- data$biomass/data$stable_AUB

ggplot(data)+
    geom_line(aes(x=year, y=AUB), color="red", size=1.5)+
    geom_line(aes(x=year, y=biomass), color="black", size=1.5)+
    geom_line(aes(x=year, y=stable_AUB), color="green", size=1.5)+
    geom_line(aes(x=year, y=dep*500000), color="red", linetype="dashed")+
    geom_line(aes(x=year, y=stable_dep*500000), color="green", linetype="dashed")+
    scale_y_continuous(limits = c(0, 550000), sec.axis = sec_axis(trans=~.*1/550000))


sim.recdevs <- read_csv(file.path(here::here(), "data_outputs", "simulation_recdevs_1.csv")) %>% summarise(across(everything(), median)) %>% select(-c(1))

pars$annual_age0devs

devs <- c(unlist(c(pars$annual_age0devs, sim.recdevs %>% as.vector)))

AUBs <- sapply(c(11:length(devs)), function(t) exp(pars$log_MeanAge0+mean(devs[(t-10):t]))*spawn.biomass.per.recruit(w=waa, m=maa, v=vaa))
plot(1989+(1:(length(devs)-10)), AUBs, type="l", ylab="Unfished Spawning Biomass (mt)", ylim=c(0, 500000), xlab="year", main="Unfished Spawning Biomass vs Fishing Mortality")

sim.recdevs <- read_csv(file.path(here::here(), "data_outputs", "simulation_recdevs.csv")) %>% select(-c(1)) %>%
                `colnames<-`(1:ncol(.)) %>%
                pivot_longer(everything(), names_to="year", values_to="rec.dev") %>%
                mutate(year = as.numeric(year))
sim.recdevs$rolling.recdev <- NA
for(i in seq(1, nrow(sim.recdevs), by=500)){
    lag <- 1
    year <- sim.recdevs[i,]$year
    for(l in 1:max(sim.recdevs$year)){
        sim.recdevs[i+l-1,]$rolling.recdev <- pars$log_MeanAge0+mean(sim.recdevs[i:(i+l-1),]$rec.dev)
    }
}

sim.aub <- sim.recdevs %>% 
                mutate(
                    AUB = exp(rolling.recdev)*spawn.biomass.per.recruit(w=waa, m=maa, v=vaa),
                    year = as.numeric(year)
                ) %>%
                group_by(year) %>%
                median_qi(AUB, .width=c(0.5, 0.95))

ggplot(sim.aub, aes(x=year, y=AUB, ymin=.lower, ymax=.upper))+
    geom_lineribbon(size=0.5)+
    scale_y_continuous(limits=c(0, 500000))+
    scale_fill_brewer(palette="Blues")
