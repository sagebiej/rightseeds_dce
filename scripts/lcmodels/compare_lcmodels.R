
library(tidyr)
library(ggplot2)
library(dplyr)

lcmodels=list()

fitstats <- data.frame(no_classes= integer() ,N= integer(), no_pars = integer(), LL = numeric(), BIC = numeric() , AIC=numeric() , AICc=numeric())

texoutput <- list()

for (i in 2:6){

  t <-  readRDS(paste0("modeloutput/lclogit/lc_",2,"class_model.rds")  )  
  
  texoutput[[i-1]] <- quicktexregapollo(t)
  
  n= t$nObs
  k= length(t$fixed)
  ll=t$LLout[1]

BIC <- -2*ll+k*log(n)
AIC <-   -2*ll + 2*k
AICc <- AIC +(2*k^2+2*k)/(n-k-1)

   
lcmodels[[paste0(i," Classes")]] <- t
  
fitstats[i-1,] <- c(i,n,k,ll , BIC, AIC , AICc)

  
}


gof_lc_graph <-fitstats %>% select(-c(N,no_pars, LL )) %>% 
  gather(key = "variable", value = "value", -no_classes) %>% 
ggplot(aes(x=no_classes, y=value)) + 
  geom_line(aes(color = variable, linetype = variable)) 






gof_lc_graph

