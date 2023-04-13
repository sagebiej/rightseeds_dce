## Estimate Mixed logit models

mixwtppath <- "modeloutput/WTPSpace/mixlogit_wtpspace/"
mixprefpath <-"modeloutput/mixlogit_prefspace/"

dir.create("modeloutput/Individual_wtp")

## change sign of cost
 
# database2 <- database2 %>% 
#   mutate(across(ends_with("x1"), ~abs(.x)*(-1)))

##### All respondents ####

database<- database2

apollo_initialise()


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_allsamples",
  modelDescr ="Simple mixedlogit model with all respondents",
  indivID    ="RID", 
 nCores=11,
  mixing = TRUE, # Mixed Logit
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody.R")


mixlog_allsamples = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                               hessianRoutine="analytic"))


apollo_saveOutput(mixlog_allsamples,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))  

#### Split1 No Text####
apollo_initialise()


database<-database2[which(database2$info_text==1),]


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_split_notext",
  modelDescr ="Simple mixedlogit model split samples No Text",
  indivID    ="RID", 
  nCores=11,
  mixing = TRUE, # Mixed Logit
  HB=FALSE,
  outputDirectory = mixwtppath
)



source("scripts/mixlogitbody.R")

mixlog_split_notext = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                               hessianRoutine="analytic"))

apollo_saveOutput(mixlog_split_notext,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))
 
individualWTP_notext <-apollo_conditionals(mixlog_split_notext, apollo_probabilities, apollo_inputs) #derive individual WTP values 

saveRDS(individualWTP_notext, file = "modeloutput/Individual_wtp/individualWTP_notext.rds")



#### Split2 Environment ####



database<-database2[which(database2$info_text==2),]


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_split_env",
  modelDescr ="Simple mixedlogit model split samples Environment",
  indivID    ="RID",
  nCores=11,
  mixing = TRUE, # Mixed Logit
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody.R")

mixlog_split_env = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))




apollo_saveOutput(mixlog_split_env,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))

individualWTP_env <-apollo_conditionals(mixlog_split_env, apollo_probabilities, apollo_inputs) #derive individual WTP values 

saveRDS(individualWTP_env, file = "modeloutput/Individual_wtp/individualWTP_env.rds")

##### Split 3 ####





database<-database2[which(database2$info_text==3),]


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_split_ind",
  modelDescr ="Simple mixedlogit model pilot split samples Industry",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody.R")

mixlog_split_ind = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))




apollo_saveOutput(mixlog_split_ind,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))


individualWTP_ind <-apollo_conditionals(mixlog_split_ind , apollo_probabilities, apollo_inputs) #derive individual WTP values 

saveRDS(individualWTP_ind, file = "modeloutput/Individual_wtp/individualWTP_ind.rds")

##### Revisions #####

### Split Sample summer/winter ###

### Estimate model for Summer


database<-database2[which(database2$jahreszeit_string=="Sommer"),]

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_season_summer",
  modelDescr ="Simple mixedlogit model split summer",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody_season.R")

mixlog_season_summer = apollo_estimate(apollo_beta, apollo_fixed,
                                   apollo_probabilities, apollo_inputs, 
                                   estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                          hessianRoutine="analytic"))




apollo_saveOutput(mixlog_season_summer,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))

### Estimate model for Winter 


database<-database2[which(database2$jahreszeit_string=="Winter"),]

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_season_winter",
  modelDescr ="Simple mixedlogit model split winter",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody_season.R")

mixlog_season_winter = apollo_estimate(apollo_beta, apollo_fixed,
                                 apollo_probabilities, apollo_inputs, 
                                 estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                        hessianRoutine="analytic"))




apollo_saveOutput(mixlog_season_winter,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))

### Split samples for different tomato types ###

### Cocktail- / Cherrytomatos / Small tomatos 

database<-database2[which(database2$q3_new==1),]

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_small_tomato",
  modelDescr ="Simple mixedlogit model small tomato",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody.R")

mixlog_small_tomato = apollo_estimate(apollo_beta, apollo_fixed,
                                       apollo_probabilities, apollo_inputs, 
                                       estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                              hessianRoutine="analytic"))




apollo_saveOutput(mixlog_small_tomato,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))

### Big tomato's 

database<-database2[which(database2$q3_new==2),]

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_big_tomato",
  modelDescr ="Simple mixedlogit model big tomato",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody.R")

mixlog_big_tomato = apollo_estimate(apollo_beta, apollo_fixed,
                                      apollo_probabilities, apollo_inputs, 
                                      estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                             hessianRoutine="analytic"))




apollo_saveOutput(mixlog_big_tomato,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))


##### Estimate standard model in preference space #####

database<-database2

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_pref_space",
  modelDescr ="Simple mixedlogit model preference space",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixprefpath
)


source("scripts/mixlogitbody_pref_space.R")

mixlog_pref_space = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))




apollo_saveOutput(mixlog_pref_space,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))


##### Estimate standard model in WTP space with fixed price


database<-database2

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_fix_price",
  modelDescr ="Simple mixedlogit model fix price",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixprefpath
)


source("scripts/mixlogitbody_fix_price.R")

mixlog_fix_price = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))




apollo_saveOutput(mixlog_fix_price,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))

##### Estimate mxl model with interaction with table 4

database<-database2 %>% filter(q14 < 4 & q15 < 4 & q16 < 4) #exlude don't knows 

modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="mixlog_int_tab4",
  modelDescr ="Simple mixedlogit model interactions table 4",
  indivID    ="RID",  mixing = TRUE, # Mixed Logit
  nCores=11,
  HB=FALSE,
  outputDirectory = mixwtppath
)


source("scripts/mixlogitbody_interaction_tab4.R")

mixlog_int_tab4 = apollo_estimate(apollo_beta, apollo_fixed,
                                   apollo_probabilities, apollo_inputs, 
                                   estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                          hessianRoutine="analytic"))




apollo_saveOutput(mixlog_int_tab4,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))