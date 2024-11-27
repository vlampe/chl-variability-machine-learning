weighted_moving_average <- function(values, weights, width){
  
  # calculates weighted moving average from vectors of
    # values: values to average
    # weights: weights ( in the form of number of observations, area... so absolute values)
  # width: window width, must be numeric and odd 
  # outut is a numeric vector
  # can handle NAs 
  
  out <- as.numeric(rep(NA, length = length(values)))
  
  lr <- (width-1)/2
  
  for(i in (lr+1):(length(values)-lr)){
    weights_i <- weights[(i-lr):(i+lr)] / sum(weights[(i-lr):(i+lr)], na.rm = T)
    
    vals_i <-  values[(i-lr):(i+lr)]
    
    mean_i <- sum((weights_i * vals_i), na.rm = T)
    #but: if all weights are 0 (no observations), make mean_i NA (not 0)
    if(all(is.na(weights_i))) mean_i <- NA
    
    out[i] <- mean_i
  }
  return(out)
}