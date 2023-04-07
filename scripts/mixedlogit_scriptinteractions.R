## Estimate Mixed logit models

No_draws =1200


## change sign of cost

# database2 <- database2 %>% 
#   mutate(across(ends_with("x1"), ~abs(.x)*(-1)))

##### All respondents ####

database<- database2

apollo_initialise()


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_scriptint",
  modelDescr ="Simple mixedlogit model with script interactions",
  indivID    ="RID", 
  nCores=11,
  mixing = TRUE, # Mixed Logit
  HB=FALSE,
  outputDirectory = "modeloutput/"
)




#rm(apollo_inputs)

randnames  <- c( "asc", "EUBIO" , "BIOVER" , "REG", "DE", "ESP" , "EIGENT", "cost"  )  
fixednames  <- c( "ascxInfo2", "EUBIOxInfo2" , "BIOVERxInfo2" , "REGxInfo2", "DExInfo2", "ESPxInfo2" , "EIGENTxInfo2", "ascxInfo3", "EUBIOxInfo3" , "BIOVERxInfo3" , "REGxInfo3", "DExInfo3", "ESPxInfo3" , "EIGENTxInfo3" )  

beta <- vector()

for (n in randnames) {
  beta[paste0("mean_",n)]=0.5
  beta[paste0("sd_",n)]=1
}

for (n in fixednames) {
  beta[n]= 0     
}


apollo_beta <-beta

# apollo_beta=c(mean_asc = 0,
#               ascwinter =-0.04,
#               mean_EUBIO = 0,
#               mean_BIOVER = 0,
#               mean_REG = 0,
#               mean_DE = 0,
#               mean_ESP = 0,
#               mean_EIGENT = 0,
#               EIGENTxOC = -0,
#               sd_asc = 1 ,
#               sd_EUBIO = 1,
#               sd_BIOVER = 1,
#               sd_REG=1,
#               sd_DE = 0.1,
#               sd_ESP = 1,
#               sd_EIGENT = 1,
#               mean_cost = 0.1,
#               sd_cost=0.4
# )

# apollo_beta=c(mean_asc = 1.2,
#               # ascwinter =-0.04,
#               # ascfemale = 0,
#               # ascage=0,
#               # aschigheredu = 0,
#               # ascincome=0,
#               mean_EUBIO = 1,
#               mean_BIOVER = 1.3,
#               mean_REG = 2.5,
#               mean_DE = 2,
#               mean_ESP = 0.5,
#               mean_EIGENT = 0.5,
#               #EIGENTxOC = -0.06,
#               sd_asc = 2 ,
#               sd_EUBIO = 1,
#               sd_BIOVER = 1,
#               sd_REG=1,
#               sd_DE = 0.5,
#               sd_ESP = 0.5,
#               sd_EIGENT = 0.4,
#               mean_cost = 0.1,
#               sd_cost=0.3,
#               mean_EUBIOX = 1,
#               mean_BIOVER = 1.3,
#               mean_REG = 2.5,
#               mean_DE = 2,
#               mean_ESP = 0.5,
#               mean_EIGENT = 0.5,
#               
# )
draws  <- paste0("draws_",randnames)
### keine Parameter fix halten
apollo_fixed <- c()

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws <- list(
  interDrawsType = "sobol",
  interNDraws    = No_draws, #50 für Code testen, min. 500 für verlässliche Ergebnisse
  interUnifDraws = c(),
  # nur normal distributed: (für cost  lognormal transformieren)
  interNormDraws = draws,
  # keine Intra-Individuen Heterogenität: (das wären abweichende Präferenzen selbe Individuen zwischen verschiedenen Choices)
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  # for (normpar in randnames) {
  # 
  #   randcoeff[[normpar]] = get(paste0("mean_",normpar)) + get(paste0("sd_",normpar)) * get(paste0("draws_",normpar))
  # 
  # }
  
  randcoeff[["asc"]] = mean_asc + sd_asc * draws_asc
  randcoeff[["EUBIO"]] = mean_EUBIO + sd_EUBIO * draws_EUBIO
  randcoeff[["BIOVER"]] = mean_BIOVER + sd_BIOVER * draws_BIOVER
  randcoeff[["REG"]] = mean_REG + sd_REG * draws_REG
  randcoeff[["DE"]] = mean_DE + sd_DE * draws_DE
  randcoeff[["ESP"]] = mean_ESP + sd_ESP * draws_ESP
  randcoeff[["EIGENT"]] = mean_EIGENT + sd_EIGENT * draws_EIGENT
  randcoeff[["cost"]] = - exp(mean_cost + sd_cost * draws_cost)
  
  
  return(randcoeff)
  
  
}

### validieren
apollo_inputs = apollo_validateInputs()


apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities (later integrated in mnl_settings below)
  V = list()
  V[['alt1']] =  -cost*(-a1_x1 +
                          asc +  ascxInfo2*(info_text==2)+ascxInfo3*(info_text==3) +
                          EUBIO*a1_EUBIO + BIOVER*a1_BIOVER + REG*a1_REG + DE*a1_DE + ESP*a1_ESP  + EIGENT*a1_x4 + 
                          EUBIOxInfo2*a1_EUBIO*(info_text==2) + BIOVERxInfo2*a1_BIOVER*(info_text==2) + REGxInfo2*a1_REG*(info_text==2) + DExInfo2*a1_DE*(info_text==2) + ESPxInfo2*a1_ESP*(info_text==2)  + EIGENTxInfo2*a1_x4*(info_text==2) +
                          EUBIOxInfo3*a1_EUBIO*(info_text==3) + BIOVERxInfo3*a1_BIOVER*(info_text==3) + REGxInfo3*a1_REG*(info_text==3) + DExInfo3*a1_DE*(info_text==3) + ESPxInfo3*a1_ESP*(info_text==3)  + EIGENTxInfo3*a1_x4*(info_text==3) 
  )
  V[['alt2']] = -cost*(-a2_x1 + 
                         asc+ ascxInfo2*(info_text==2)+ascxInfo3*(info_text==3) +
                         EUBIO*a2_EUBIO + BIOVER*a2_BIOVER  + REG*a2_REG + DE*a2_DE  + ESP*a2_ESP  + EIGENT*a2_x4 +
                         EUBIOxInfo2*a2_EUBIO*(info_text==2) + BIOVERxInfo2*a2_BIOVER*(info_text==2) + REGxInfo2*a2_REG*(info_text==2) + DExInfo2*a2_DE*(info_text==2) + ESPxInfo2*a2_ESP*(info_text==2)  + EIGENTxInfo2*a2_x4*(info_text==2) +
                         EUBIOxInfo3*a2_EUBIO*(info_text==3) + BIOVERxInfo3*a2_BIOVER*(info_text==3) + REGxInfo3*a2_REG*(info_text==3) + DExInfo3*a2_DE*(info_text==3) + ESPxInfo3*a2_ESP*(info_text==3)  + EIGENTxInfo3*a2_x4*(info_text==3) 
  
                       )  
  V[['alt3']] = 0
  
  
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3) ,
    avail         = 1, # all alternatives are available in every choice
    choiceVar     = pref1,
    V             = V  # tell function to use list vector defined above
    
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws - nur bei Mixed Logit!
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}




mixlog_scriptint = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))

apollo_saveOutput(mixlog_scriptint,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))
