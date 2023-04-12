apollo_beta=c(asc = 4,
              ascwinter =0,
              cost = -0.7,
              EUBIO = 1.17,
              BIOVER = 1.5,
              REG = 2.5 ,
              DE = 2,
              ESP = 0.5,
              EIGENT = 0.5,
              EIGENTxOC = 0
)


### keine Parameter fix halten
apollo_fixed = c()

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
  V[['alt1']] =   -cost*(asc -a1_x1 + ascwinter*winter + EUBIO*a1_EUBIO + BIOVER*a1_BIOVER +REG*a1_REG + DE*a1_DE + ESP*a1_ESP  + EIGENT*a1_x4+EIGENTxOC*a1_x4*oc)
  V[['alt2']] =   -cost*(asc -a2_x1 + ascwinter*winter + EUBIO*a2_EUBIO + BIOVER*a2_BIOVER + REG*a2_REG+ DE*a2_DE  + ESP*a2_ESP + EIGENT*a2_x4+EIGENTxOC*a2_x4*oc)
  V[['alt3']] = 0
  
  # V[['alt1']] =  asc +  ascwinter*winter + cost*a1_x1 + EUBIO*a1_EUBIO + BIOVER*a1_BIOVER + REG*a1_REG + DE*a1_DE + ESP*a1_ESP + EIGENT*a1_x4 + EIGENTxOC*a1_x4*oc
  # 
  # V[['alt2']] =  asc +  ascwinter*winter + cost*a2_x1  + EUBIO*a2_EUBIO + BIOVER*a2_BIOVER + REG*a2_REG + DE*a2_DE + ESP*a2_ESP + EIGENT*a2_x4 + EIGENTxOC*a2_x4*oc
  
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
  ### P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

