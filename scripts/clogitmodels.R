## Estimate clogit model

# database2 <- database2 %>% 
#   mutate(across(ends_with("x1"), ~abs(.x)*(-1)))


##### All respondents ####

database <- database2

apollo_initialise()


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="clog_allsamples",
  modelDescr ="Simple MNL model with all respondents",
  indivID    ="RID", 
  nCores=11,
  outputDirectory = "modeloutput/"
)


source("scripts/clogitbody.R")


clog_allsamples = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                               hessianRoutine="analytic"))





apollo_saveOutput(clog_allsamples,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))  

#### Split1 No Text####



database<-database2[which(database2$info_text==1),]


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="clog_split_notext",
  modelDescr ="Simple MNL model pilot split samples No Text",
  indivID    ="RID",
  outputDirectory = "modeloutput/"
)



source("scripts/clogitbody.R")

clog_split_notext = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                               hessianRoutine="analytic"))




apollo_saveOutput(clog_split_notext,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))
 




#### Split2 Environment ####



database<-database2[which(database2$info_text==2),]


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="clog_split_env",
  modelDescr ="Simple MNL model pilot split samples Environment",
  indivID    ="RID",
  outputDirectory = "modeloutput/"
)


source("scripts/clogitbody.R")

clog_split_env = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))




apollo_saveOutput(clog_split_env,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))



##### Split 3 ####





database<-database2[which(database2$info_text==3),]


modelOutput_settings = list(printPVal=T)

### Set core controls
apollo_control = list(
  modelName  ="clog_split_ind",
  modelDescr ="Simple MNL model pilot split samples Industry",
  indivID    ="RID",
  outputDirectory = "modeloutput/"
)


source("scripts/clogitbody.R")

clog_split_ind = apollo_estimate(apollo_beta, apollo_fixed,
                                    apollo_probabilities, apollo_inputs, 
                                    estimate_settings=list(maxIterations=400, estimationRoutine="bfgs",
                                                           hessianRoutine="analytic"))




apollo_saveOutput(clog_split_ind,  saveOutput_settings = list(saveCov=FALSE, saveCorr=FALSE))




