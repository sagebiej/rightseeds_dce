

rm(list = ls(, pattern = "^apollo_"))





database <- database2



fixednames  <- c("asc"  , "cost" , "EUBIO" , "BIOVER"  ,"REG" ,  "DE",    "ESP",   "EIGENT" ,    "EIGENTxOC","delta" )


apollo_initialise()

### Set core controls
apollo_control = list(
  modelName  =paste0("lc_",n_classes,"class"),
  modelDescr =paste0(n_classes," Class LCL model for Tomato DCE"),
  indivID    ="RID",
  nCores     = 1, 
    outputDirectory = "modeloutput/lclogit/"
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #
### Vector of parameters, including any that are kept fixed in estimation

beta <-as.vector(t(outer(fixednames, letters[seq_along(1:n_classes)], paste, sep="_"))) 
apollo_beta= rep(0,length(beta))
names(apollo_beta) <-beta

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("delta_a")

# ################################################################# #
#### DEFINE LATENT CLASS COMPONENTS                              ####
# ####################################################  
 
apollo_lcPars=function(apollo_beta, apollo_inputs){


 attrib  <- fixednames[-grep("delta",fixednames)]
  lcpars <- list()
  
  for (f in fixednames) {
    lcpars[[f]] <- map(seq_along(1:n_classes), ~ get(paste0(f, "_", letters[.])))
    
  }

  V=list()
  
  for (c in letters[seq_along(1:n_classes)]) {
    V[paste0("class_",c)] = get(paste0("delta_",c))
    
  }

  
  mnl_settings = list(
    alternatives = set_names(1:n_classes, paste0("class_",letters[1:n_classes])), 
    avail        = 1, 
    choiceVar    = NA, 
    V            = V
  )
  
  lcpars[["pi_values"]] = apollo_mnl(mnl_settings, functionality="raw")
    lcpars[["pi_values"]] = apollo_firstRow(lcpars[["pi_values"]], apollo_inputs)
  
  return(lcpars)
 

}


# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  
  
  ### Define settings for MNL model component that are generic across classes
  mnl_settings = list(
    alternatives = c(alt1=1, alt2=2, alt3=3),
    avail        = 1,
    choiceVar    = pref1
  )
  
  
  ### Loop over classes
  
  
  for (s in seq_along(1:n_classes)) {
    
    ### Compute class-specific utilities
     V=list()


    # V[['alt1']] = -cost[[s]]*(-a1_x1 +
    #                        asc[[s]] +
    #                        EUBIO[[s]]*a1_EUBIO + BIOVER[[s]]*a1_BIOVER + REG[[s]]*a1_REG + DE[[s]]*a1_DE + ESP[[s]]*a1_ESP  + EIGENT[[s]]*a1_x4 +EIGENTxOC[[s]]*a1_x4*oc)
    # 
    #     V[['alt2']] = -cost[[s]]*(-a2_x1 +
    #                           asc[[s]] +
    #                           EUBIO[[s]]*a2_EUBIO + BIOVER[[s]]*a2_BIOVER + REG[[s]]*a2_REG + DE[[s]]*a2_DE + ESP[[s]]*a2_ESP  + EIGENT[[s]]*a2_x4 +EIGENTxOC[[s]]*a2_x4*oc)

     
     
     V[['alt1']] = cost[[s]]*a1_x1 +
                                 asc[[s]] +
                                 EUBIO[[s]]*a1_EUBIO + BIOVER[[s]]*a1_BIOVER + REG[[s]]*a1_REG + DE[[s]]*a1_DE + ESP[[s]]*a1_ESP  + EIGENT[[s]]*a1_x4 +EIGENTxOC[[s]]*a1_x4*oc

     V[['alt2']] = cost[[s]]*a2_x1 +
                                 asc[[s]] +
                                 EUBIO[[s]]*a2_EUBIO + BIOVER[[s]]*a2_BIOVER + REG[[s]]*a2_REG + DE[[s]]*a2_DE + ESP[[s]]*a2_ESP  + EIGENT[[s]]*a2_x4 +EIGENTxOC[[s]]*a2_x4*oc


    V[['alt3']] = 0 
    
    
    mnl_settings$V = V
    mnl_settings$componentName = paste0("Class_",s)
    
    ### Compute within-class choice probabilities using MNL model
    P[[paste0("Class_",s)]] = apollo_mnl(mnl_settings, functionality)
    
    ### Take product across observation for same individual
    P[[paste0("Class_",s)]] = apollo_panelProd(P[[paste0("Class_",s)]], apollo_inputs ,functionality)
    
    
  }
  
  ### Compute latent class model probabilities
  lc_settings   = list(inClassProb = P, classProb=pi_values)
  P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
} 


lc2cl = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(hessianRoutine="maxLik"))


lc2cl[["unconditionals"]] <- apollo_lcUnconditionals(lc2cl,apollo_probabilities,apollo_inputs)

apollo_saveOutput(lc2cl,  saveOutput_settings = list(saveEst=TRUE, saveCov=FALSE, saveCorr=FALSE))





