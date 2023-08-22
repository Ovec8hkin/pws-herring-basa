# init.admb.params <- function(reps){
    
#     pars <- r4ss::read.admbFit("PWS_ASA")
#     inits <- list()
#     for(j in 1:reps){
#         # This is a check on the CV of pars - basically high uncertainty indicates parameter is fixed or uninformed
#         # CV near 0 indicates fixed parameter (not estimated), CV near inf indicates uninformed parameter.
#         # Use CV = 0.5 for drawing parameters (arbitrary) if parameter is fixed or uninformed.
#         # This may cause boundary issues in NUTS
#         # if(abs(pars$std[1]/pars$est[1])<1 | !is.infinite(pars$std[1]/pars$est[1])){
#         #     inits[[j]] <- rnorm(1, pars$est[1], sd=pars$std[1])
#         # }else{
#         #     inits[[j]] <- rnorm(1, pars$est[1], sd=0.5 * pars$est[1])
#         # }

#         sd <- ifelse(abs(pars$std[1]/pars$est[1])<1 | !is.infinite(pars$std[1]/pars$est[1]), pars$std[1], pars$est[1]/2)

#         inits[[j]] <- rnorm(1, pars$est[1], sd=sd)
        
#         for(k in 2:length(pars$est)){
#             # Checks for parameters bounded between 0 and 1 (assumes these have "prob" in name)
#             if((grepl("prob", pars$names[k])) & pars$est[k]<=1 & pars$est[k]>=0){
#                 inits[[j]] <- c(inits[[j]], runif(1, 0.01, 0.99))
#             #inits[[j]] <- c(inits[[j]],rbeta(1,pars$est[k],0.99))
#             }else if(abs(pars$std[k]/pars$est[k])<1 | !is.infinite(pars$std[k]/pars$est[k])){
#                 inits[[j]] <- c(inits[[j]], rnorm(1, pars$est[k], sd=pars$std[k]))  
#             }else{
#                 inits[[j]] <- c(inits[[j]], rnorm(1, pars$est[k], sd=0.5 * pars$est[k]))
#             }
#         }
#         # Set names
#         names(inits[[j]]) <- pars$names
#         # Omit derived quantities included in the ADMB fit output (e.g. from sdreport)
#         rem <- which(pars$names=="SSB_final_year")
#         inits[[j]] <- inits[[j]][-rem]
#     }
#     return(inits)
# }

init.admb.params <- function(reps){
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
    return(inits)
}