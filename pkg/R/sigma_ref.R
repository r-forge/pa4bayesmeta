
# computation of the reference standard deviation
sigma_ref <- function(df, 
                      type.sigma.ref = "geometric") {
  # input:
  # df: data frame with one column "sigma" containing the standard deviations sigmai in each study
  # type.sigma.ref: type of computation of the sigma_ref:
  ## "geometric" geometric mean of within-study sds as suggested in equation (7) by Sorbye and Rue (scaling)
  ## "harmonic" weighted harmonic mean as suggested on p. 490 by Hoaglin (2016) "Misunderstandings about Q..."
  # output:
  # reference standard deviation
  
  if(! type.sigma.ref %in% c("geometric", "harmonic"))
    stop("The specified type.sigma.ref is not supported. Possible values are 'geometric' and 'harmonic'.")
  
  sigma <- df$sigma
  
  if (type.sigma.ref == "geometric") {
    val <- exp(mean(log(sigma)))
  }
  
  if (type.sigma.ref == "harmonic")
  {
    # number of studies
    kk <- length(sigma)
    # weights
    wi <- 1 / (sigma ^ 2)
    # variance
    ss2 <- ((kk - 1) * sum(wi)) / ((sum(wi)) ^ 2 - sum(wi ^ 2))
    # standard deviation
    val <- sqrt(ss2)
  }
  
  return(val)
}