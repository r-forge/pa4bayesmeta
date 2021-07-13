
pri_par_adjust <- function(df,
                           rlmc = 0.5,
                           tail.prob = 0.5,
                           distributions = c("HN", "HC"),
                           type.sigma.ref = "geometric") {
  # function for a dynamic median RLMC-based scaling parameters adjustment for HN, EXP, HC, LMX
  
  # check parameter values
  if(! (rlmc >= 0 & rlmc <= 1))
    stop("rlmc must be in [0, 1].")
  if(! (tail.prob >= 0 & tail.prob <= 1))
    stop("tail.prob must be in [0, 1].")
  if(! type.sigma.ref %in% c("geometric", "harmonic"))
    stop("The specified type.sigma.ref is not supported. Possible values are 'geometric' and 'harmonic'.")
  
  # supporting functions
  
  # ### Analytical formulae to find the scaling factor for a prior distribution fulfilling tail-adjustment given a threshold UU and an alpha tail-probability
  # ### U and alpha fixed -> how much is AA?
  
  AA_from_Ualpha_HN <- function(UU, alpha) {
    return(UU / qnorm(
      1 - alpha / 2,
      mean = 0,
      sd = 1,
      lower.tail = TRUE
    ))
  }
  
  AA_from_Ualpha_Exp <- function(UU, alpha) {
    return(-UU / log(alpha))
  }
  
  AA_from_Ualpha_HC <- function(UU, alpha) {
    return(UU / tan(pi * (1 - alpha) / 2))
  }
  
  AA_from_Ualpha_Lomax <- function(UU, alpha) {
    return(UU * alpha / (1 - alpha))
  }
  
  
  # computation of the reference threshold U_ref
  # P[tau>U_ref]=tail.prob
  U_ref <-
    sqrt(rlmc / (1 - rlmc)) * sigma_ref(df = df, type.sigma.ref = type.sigma.ref)
  
  parameters <- list()
  # tail.prob adjusting
  for (i in seq_along(distributions)) {
    if (!distributions[i] %in% c("HN", "HC", "EXP", "LMX"))
      stop(paste0(
        "The distribution ",
        distributions[i],
        " is not supported by this function."
      ))
    
    if (distributions[i] == "HN") {
      p_HN <- AA_from_Ualpha_HN(U_ref, alpha = tail.prob)
      parameters[[i]] <- p_HN
      names(parameters)[i] <- "p_HN"
    }
    if (distributions[i] == "HC") {
      p_HC <- AA_from_Ualpha_HC(U_ref, alpha = tail.prob)
      parameters[[i]] <- p_HC
      names(parameters)[i] <- "p_HC"
    }
    if (distributions[i] == "EXP") {
      p_EXP <- AA_from_Ualpha_Exp(U_ref, alpha = tail.prob)
      parameters[[i]] <- p_EXP
      names(parameters)[i] <- "p_EXP"
    }
    if (distributions[i] == "LMX") {
      p_LMX <- AA_from_Ualpha_Lomax(U_ref, alpha = tail.prob)
      parameters[[i]] <- p_LMX
      names(parameters)[i] <- "p_LMX"
    }
  }
  
  return(parameters)
}
