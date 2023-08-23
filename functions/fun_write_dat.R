# fun_write_dat.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Project population dynamics to next year with given catches

fun_write_dat <- function(dat.files, years.to.remove){
  
  options(scipen=10)

  # Write main data file with actual catches (known w/ error)
  # and observed survey values (from fun_obsm.R).
  PWS.ASA.dat <- dat.files$PWS_ASA.dat
  # Add new data 
  PWS.ASA.dat[[1]]  <- PWS.ASA.dat[[1]]-years.to.remove                                              # Number of year in data
  PWS.ASA.dat[[2]]  <- PWS.ASA.dat[[2]]-years.to.remove                                              # Number of years to be fit to
  # 3 is number of age classes (10). Needs no change.
  PWS.ASA.dat[[4]]  <- PWS.ASA.dat[[4]][1:(nrow(PWS.ASA.dat[[4]])-years.to.remove),]                 # Weight-at-age 
  PWS.ASA.dat[[5]]  <- PWS.ASA.dat[[5]][1:(nrow(PWS.ASA.dat[[5]])-years.to.remove),]                 # Fecundity
  PWS.ASA.dat[[6]]  <- PWS.ASA.dat[[6]][1:(nrow(PWS.ASA.dat[[6]])-years.to.remove),]                 # Pound catch fishery
  # 7 is prop of pound-catch killed (0.25). Needs no change.
  PWS.ASA.dat[[8]]  <- PWS.ASA.dat[[8]][1:(nrow(PWS.ASA.dat[[8]])-years.to.remove),]                 # Foodbait fishery
  PWS.ASA.dat[[9]]  <- PWS.ASA.dat[[9]][1:(nrow(PWS.ASA.dat[[9]])-years.to.remove),]                 # Gillnet fishery
  PWS.ASA.dat[[10]] <- as.matrix(PWS.ASA.dat[[10]][1:(nrow(PWS.ASA.dat[[10]])-years.to.remove),])    # Seine net fishery
  PWS.ASA.dat[[11]] <- as.matrix(PWS.ASA.dat[[11]][1:(nrow(PWS.ASA.dat[[11]])-years.to.remove),])    # % female

  PWS.ASA.dat[[12]] <- as.matrix(PWS.ASA.dat[[12]][1:(nrow(PWS.ASA.dat[[12]])-years.to.remove),])    # MDM
  PWS.ASA.dat[[13]] <- as.matrix(PWS.ASA.dat[[13]][1:(nrow(PWS.ASA.dat[[13]])-years.to.remove),])    # egg deposition
  PWS.ASA.dat[[14]] <- as.matrix(PWS.ASA.dat[[14]][1:(nrow(PWS.ASA.dat[[14]])-years.to.remove),])    # egg deposition sd
  # 15 is first year of ADFG hydro survey. Needs no change.
  PWS.ASA.dat[[16]] <- as.matrix(PWS.ASA.dat[[16]][1:(nrow(PWS.ASA.dat[[16]])-years.to.remove),])    # ADFG hydroacoustic survey biomass
  # 17 is first year of PWSSC hydro survey. Needs no change.
  PWS.ASA.dat[[18]] <- as.matrix(PWS.ASA.dat[[18]][1:(nrow(PWS.ASA.dat[[18]])-years.to.remove),])    # PWSSC hydroacoustic survey biomass
  PWS.ASA.dat[[19]] <- as.matrix(PWS.ASA.dat[[19]][1:(nrow(PWS.ASA.dat[[19]])-years.to.remove),])    # PWSSC hydroacoustic survey biomass sd
  PWS.ASA.dat[[20]] <- PWS.ASA.dat[[20]][1:(nrow(PWS.ASA.dat[[20]])-years.to.remove),]               # seine age composition
  PWS.ASA.dat[[21]] <- PWS.ASA.dat[[21]][1:(nrow(PWS.ASA.dat[[21]])-years.to.remove),]               # spawner age composition
  PWS.ASA.dat[[22]] <- as.matrix(PWS.ASA.dat[[22]][1:(nrow(PWS.ASA.dat[[22]])-years.to.remove),])    # juvenile aerial survey
  PWS.ASA.dat[[24]] <- PWS.ASA.dat[[24]][1:(nrow(PWS.ASA.dat[[24]])-years.to.remove),]               # spawner age composition
  

  write.data(PWS.ASA.dat, "PWS_ASA.dat")

  # Write ESS data file with actual sample sizes (fixed).
  PWS.ASA.ESS.ctl <- dat.files$PWS_ASA_ESS.ctl 
  PWS.ASA.ESS.ctl[[1]] <- as.matrix(PWS.ASA.ESS.ctl[[1]][1:(nrow(PWS.ASA.ESS.ctl[[1]])-years.to.remove),])
  PWS.ASA.ESS.ctl[[2]] <- as.matrix(PWS.ASA.ESS.ctl[[2]][1:(nrow(PWS.ASA.ESS.ctl[[2]])-years.to.remove),])
  PWS.ASA.ESS.ctl[[3]] <- as.matrix(PWS.ASA.ESS.ctl[[3]][1:(nrow(PWS.ASA.ESS.ctl[[3]])-years.to.remove),])
  # Write file without headers
  write.data(PWS.ASA.ESS.ctl, "PWS_ASA(ESS).ctl")
  
  
  # Write covariate file (not currently used)
  PWS.ASA.covariate.ctl <- dat.files$PWS_ASA_covariate.ctl
  # Skip indices 1-4 for now, as they do not change with year
  PWS.ASA.covariate.ctl[[5]] <- as.matrix(PWS.ASA.covariate.ctl[[5]][1:(nrow(PWS.ASA.covariate.ctl[[5]])-years.to.remove),]) 
  PWS.ASA.covariate.ctl[[6]] <- as.matrix(PWS.ASA.covariate.ctl[[6]][1:(nrow(PWS.ASA.covariate.ctl[[6]])-years.to.remove),])
  # Skip indices 7-11 for now, as they do not change with year
  PWS.ASA.covariate.ctl[[12]] <- PWS.ASA.covariate.ctl[[12]][1:(nrow(PWS.ASA.covariate.ctl[[12]])-years.to.remove),]
  PWS.ASA.covariate.ctl[[14]] <- as.matrix(PWS.ASA.covariate.ctl[[14]][1:(nrow(PWS.ASA.covariate.ctl[[14]])-years.to.remove),])
  
  # Write .dat file without headers
  write.data(PWS.ASA.covariate.ctl, "PWS_ASA(covariate).ctl")
  
  # Write age comp sample sizes file with actual 
  # sample sizes (currently fixed).
  age.comp.samp.sizes.txt <- dat.files$agecomp_samp_sizes.txt
  age.comp.samp.sizes.txt[[1]] <- as.matrix(age.comp.samp.sizes.txt[[1]][1:(nrow(age.comp.samp.sizes.txt[[1]])-years.to.remove),])
  age.comp.samp.sizes.txt[[2]] <- as.matrix(age.comp.samp.sizes.txt[[2]][1:(nrow(age.comp.samp.sizes.txt[[2]])-years.to.remove),])
  age.comp.samp.sizes.txt[[3]] <- as.matrix(age.comp.samp.sizes.txt[[3]][1:(nrow(age.comp.samp.sizes.txt[[3]])-years.to.remove),])
  
  # Write file without headers
  write.table(
      rbind(" ", " ", " ", " ", age.comp.samp.sizes.txt[[1]], " ", " "),
      file = "agecomp_samp_sizes.txt", 
      append = F, sep = " ",
      row.names=FALSE, col.names=FALSE,
      quote=F
    )
    for(i in 2:length(age.comp.samp.sizes.txt)){
        write.table(
          rbind(age.comp.samp.sizes.txt[[i]], " ", " "),
          file = "agecomp_samp_sizes.txt", 
          append = T, sep = " ",
          row.names=FALSE, col.names=FALSE,
          quote=F
        )
    }
  
  PWS.ASA.disease.dat <- dat.files$PWS_ASA_disease.dat
  PWS.ASA.disease.dat[[1]] <- PWS.ASA.disease.dat[[1]][1:(nrow(PWS.ASA.disease.dat[[1]])-years.to.remove),]  # Observed VHSV antibodies
  # Skip indices 2-4, as they do not change with year.
  PWS.ASA.disease.dat[[5]] <- PWS.ASA.disease.dat[[5]][1:(nrow(PWS.ASA.disease.dat[[5]])-years.to.remove),]   # Observed ichthyophonus antibodies
  # Skip indices 6-8, as they do not change with year.
  write.data(PWS.ASA.disease.dat, "PWS_ASA_disease.dat")

  options(scipen=0)

}

write.data <- function(data, fname){
    write.table(
      rbind(" ", data[[1]], " "),
      file = fname, 
      append = F, sep = " ",
      row.names=FALSE, col.names=FALSE,
      quote=F
    )
    for(i in 2:length(data)){
        write.table(
          rbind(data[[i]], " "),
          file = fname, 
          append = T, sep = " ",
          row.names=FALSE, col.names=FALSE,
          quote=F
        )
    }
}