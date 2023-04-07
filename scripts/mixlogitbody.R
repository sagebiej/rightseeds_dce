

#rm(apollo_inputs)

randnames  <- c( "asc", "EUBIO" , "BIOVER" , "REG", "DE", "ESP" , "EIGENT", "cost"  )  
fixednames  <- c( "ascwinter" , "EIGENTxOC")  

beta <- vector()

for (n in randnames) {
  beta[paste0("mean_",n)]=0.5
  beta[paste0("sd_",n)]=1
}

for (n in fixednames) {
  beta[n]= 0     
}

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

apollo_beta=c(mean_asc = 1.2,
              ascwinter =-0.04,
              ascfemale = 0,
              ascage=0,
              aschigheredu = 0,
              ascincome=0,
              mean_EUBIO = 1,
              mean_BIOVER = 1.3,
              mean_REG = 2.5,
              mean_DE = 2,
              mean_ESP = 0.5,
              mean_EIGENT = 0.5,
              EIGENTxOC = -0.06,
              sd_asc = 2 ,
              sd_EUBIO = 1,
              sd_BIOVER = 1,
              sd_REG=1,
              sd_DE = 0.5,
              sd_ESP = 0.5,
              sd_EIGENT = 0.4,
              mean_cost = 0.1,
              sd_cost=0.3
)
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
                          asc + ascwinter*winter + ascfemale*Female_s + ascage*Age_s + aschigheredu*HigherEdu_s + ascincome*Income_s + 
                          EUBIO*a1_EUBIO + BIOVER*a1_BIOVER + REG*a1_REG + DE*a1_DE + ESP*a1_ESP  + EIGENT*a1_x4 +EIGENTxOC*a1_x4*oc)
  V[['alt2']] = -cost*(-a2_x1 + 
                        asc+ ascwinter*winter + ascfemale*Female_s + ascage*Age_s + aschigheredu*HigherEdu_s + ascincome*Income_s +
                         EUBIO*a2_EUBIO + BIOVER*a2_BIOVER  + REG*a2_REG + DE*a2_DE  + ESP*a2_ESP  + EIGENT*a2_x4 +EIGENTxOC*a2_x4*oc)  
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

