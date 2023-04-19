
calculate_avg_lc_wtp <- function(lcmodels, lcmodelsWTP, num_classes) {
  # Extract probability values
  prob <- numeric(length = num_classes)
  for (i in 1:num_classes) {
    prob[i] <- as.numeric(str_extract(lcmodels$componentReport$model$para[[3+i]], "\\d+\\.\\d+"))
  }
  
  # Calculate WTP for all attributes
  num_attributes <- length(lcmodelsWTP[[1]]$wtp)
  wtp <- numeric(length = num_attributes)
  for (j in 1:num_attributes) {
    wtp_att <- 0
    for (i in 1:num_classes) {
      wtp_att <- wtp_att + prob[i] * lcmodelsWTP[[paste0("Class ", i)]]$wtp[[j]]
    }
    wtp[j] <- wtp_att
  }
  
  return(wtp)
}






