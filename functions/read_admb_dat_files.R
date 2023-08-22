read.admb.files <- function(){
    # Parameters that are fixed and thus included as data
    #Z_3_8 <- 0.25
    egg.add <- 0.4
    
    PWS.ASA.dat <- data.reader(filename="PWS_ASA.dat") # This is nyr - we want to start at nyr_tobefit
    nyr.tot           <- PWS.ASA.dat[[1]]
    PWS.ASA.dat       <- PWS.ASA.dat[-1]
    nyr.fit           <- PWS.ASA.dat[[1]]
    nage              <- PWS.ASA.dat[[2]]
    weight.at.age     <- PWS.ASA.dat[[3]] 
    fecundity         <- PWS.ASA.dat[[4]]
    pound.catch       <- PWS.ASA.dat[[5]]
    pk                <- PWS.ASA.dat[[6]]
    food.bait.catch   <- PWS.ASA.dat[[7]]
    gillnet.catch     <- PWS.ASA.dat[[8]]
    seine.catch       <- PWS.ASA.dat[[9]]
    f.female.spawners <- PWS.ASA.dat[[10]]
    
    mdm.survey            <- PWS.ASA.dat[[11]]
    eggdep.survey         <- PWS.ASA.dat[[12]]
    eggdep.cv             <- PWS.ASA.dat[[13]]
    adfg.hydro.survey     <- PWS.ASA.dat[[15]]
    pwssc.hydro.survey    <- PWS.ASA.dat[[17]]
    pwssc.hydro.cv        <- PWS.ASA.dat[[18]]
    seine.age.comp        <- PWS.ASA.dat[[19]]
    spawn.age.comp        <- PWS.ASA.dat[[20]]
    juvenile.survey       <- PWS.ASA.dat[[21]]
    sero.comp             <- PWS.ASA.dat[[23]]
    
    PWS.ASA.ESS.ctl <- data.reader(filename="PWS_ASA(ESS).ctl")
    seine.ess <- PWS.ASA.ESS.ctl[[1]]
    spawn.ess <- PWS.ASA.ESS.ctl[[2]]
    sero.ess  <- PWS.ASA.ESS.ctl[[3]]
    
    seine.indices         <- which(rowSums(seine.age.comp[1:nyr.tot, ])>0)
    spawnsurvey.indices   <- which(rowSums(spawn.age.comp[1:nyr.tot, ])>0)
    egg.indices           <- which(eggdep.survey[1:nyr.tot]>0)
    hydADFG.indices       <- which(adfg.hydro.survey[1:nyr.tot]>0)
    hydPWSSC.indices      <- which(pwssc.hydro.survey[1:nyr.tot]>0)
    juvenile.indices      <- which(juvenile.survey[1:nyr.tot]>0)
    sero.indices           <- which(apply(sero.comp, 2, function(x) any(x>=0)))
    
    model.data <- list(nyr.tot=nyr.tot,
                        nyr.fit=nyr.fit,
                        nage=nage,
                        weight.at.age=weight.at.age[1:nyr.tot, ],
                        fecundity=fecundity[1:nyr.tot, ],
                        pound.catch=pound.catch[1:nyr.tot, ],
                        pk=pk,
                        food.bait.catch=food.bait.catch[1:nyr.tot, ],
                        gillnet.catch=gillnet.catch[1:nyr.tot, ],
                        seine.catch=seine.catch[1:nyr.tot],
                        f.female.spawners=f.female.spawners[1:nyr.tot],
                        seine.ess=seine.ess,
                        spawn.ess=spawn.ess,
                        sero.ess=sero.ess,
                        seine.age.comp=seine.age.comp[1:nyr.tot, ],
                        spawn.age.comp=spawn.age.comp[1:nyr.tot, ],
                        mdm.survey=mdm.survey[1:nyr.tot],
                        eggdep.survey=eggdep.survey[1:nyr.tot],
                        eggdep.cv=eggdep.cv[1:nyr.tot],
                        adfg.hydro.survey=adfg.hydro.survey[1:nyr.tot],
                        pwssc.hydro.survey=pwssc.hydro.survey[1:nyr.tot],
                        pwssc.hydro.cv=pwssc.hydro.cv[1:nyr.tot],
                        seine.indices=seine.indices,
                        spawnsurvey.indices=spawnsurvey.indices,
                        egg.indices=egg.indices,
                        hydADFG.indices=hydADFG.indices,
                        hydPWSSC.indices=hydPWSSC.indices,
                        #Z.3.8=Z.3.8,
                        egg.add=egg.add,
                        juvenile.survey=juvenile.survey[1:nyr.tot],
                        juvenile.indices=juvenile.indices,
                        sero.comp = sero.comp[1:nyr.tot,],
                        sero.indices = sero.indices)
    return(model.data)
}
