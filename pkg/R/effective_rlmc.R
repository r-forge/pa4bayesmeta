
####---- effective Relative Latent Model Complexity computation ----####

effective_rlmc <- function(df,
                           r.tau.prior,
                           MM = 10 ^ 6,
                           output = "sample",
                           step = ifelse(output == "prob", 0.03, NULL)) {
  # computation of the effective relative latent model complexity by MC sampling
  # input:
  # df: data frame containing a column df$sigma
  # r.tau.prior: randomisation function for the prior
  # MM: number of MC samples
  # output: "sample", "summary", "prob"
  # step: step value for bins for "prob"
  # results:
  # descriptive statistics or sample for the effective MRLMC
  # sample: sample
  # summary: descriptive statistics: summary: min, 0.25, 0.5, mean, 0.75, max
  # prob: a data frame with x,y columns for plotting
  
  if(! output %in% c("sample", "summary", "prob"))
    stop("The specified output is not supported. Possible values are 'sample', 'summary' and 'prob'.")
  
  set.seed(12567)
  
  tau_sim <- r.tau.prior(MM)
  kk <- length(df$sigma) # number od studies in the data frame
  pdsum <- 0
  for (i in 1:kk) {
    sim_ICCi <- tau_sim ^ 2 / (tau_sim ^ 2 + df$sigma[i] ^ 2)
    pdsum <- pdsum + sim_ICCi
  }
  if (output == "sample") {
    return(pdsum / kk)
  }
  if (output == "summary") {
    return(summary(pdsum / kk))
  }
  if (output == "prob") {
    data <- pdsum / kk
    breaks <- seq(0, 1, by = step)
    bin <- cut(data, breaks, include.lowest = TRUE)
    est <- tabulate(bin, length(levels(bin)))
    y <- est / (diff(breaks) * length(data))
    x <- breaks[-1] - step / 2
    return(data.frame(x = x, y = y))
  }
}