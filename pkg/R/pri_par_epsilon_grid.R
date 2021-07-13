
####---- epsilon-grid for scaled A*|X| HN, EXP, HC, LMX heterogeneity priors ----####


pri_par_epsilon_grid <- function(distributions = c("HN", "HC"),
                                 AA0.HN = 1,
                                 AA0.HC = 1,
                                 AA0.EXP = NULL,
                                 AA0.LMX = NULL,
                                 grid.epsilon = 0.00354) {
  # Function for automatic epsilon-grid computation for HN, EXP, HC, LMX applying the A*|X| methodology
  # input:
  # AA0.HN (EXP, HC, LMX): A*|X| scaling of the base heterogeneity prior HN (EXP, HC, LMX)
  # grid.epsilon: value of the grid.epsilon for epsilon-grid (only two values lower and upper) computation
  # see Roos et al. (2015) for details
  # output:
  # lower and upper parameters for HN, EXP, HC, LMX
  
  # if(! all(distributions %in% c("HN", "HC", "EXP", "LMX")))
  #   stop(paste0("At least one of the specified distributions is not supported. ",
  #        "The supported distributions are 'HN', 'HC', 'EXP' and 'LMX'."))
  
  if(! (is.null(AA0.HN) || AA0.HN > 0))
    stop("AA0.HN must either be a positive real number or NULL.")
  if(! (is.null(AA0.HC) || AA0.HC > 0))
    stop("AA0.HC must either be a positive real number or NULL.")
  if(! (is.null(AA0.EXP)|| AA0.EXP > 0))
    stop("AA0.EXP must either be a positive real number or NULL.")
  if(! (is.null(AA0.LMX)|| AA0.LMX > 0))
    stop("AA0.LMX must either be a positive real number or NULL.")
  
  if(! grid.epsilon > 0)
    stop("grid.epsilon must be a positive real number.")
  
  # supporting functions
  
  
  # HN
  HN_A0_2_Al_Au <- function(AA0, eps = grid.epsilon) {
    ## AA0: scaling of the base half normal distribution
    ## eps: epsilon for the epsilon local grid
    ## output AAl, AAu: epsilon local grid for the scaled half normal distribution
    AAl <-
      AA0 * (1 / (1 - eps ^ 2) ^ 2 - sqrt(1 / (1 - eps ^ 2) ^ 4 - 1))
    AAu <-
      AA0 * (1 / (1 - eps ^ 2) ^ 2 + sqrt(1 / (1 - eps ^ 2) ^ 4 - 1))
    return(c(AAl, AAu))
  }
  
  # EXP
  EXP_A0_2_Al_Au <- function(AA0, eps = grid.epsilon) {
    ## AA0: scaling of the base exponential distribution
    ## eps: epsilon for the epsilon local grid
    ## output AAl, AAu: epsilon local grid for the scaled exponential distribution
    AAl <-
      2 * AA0 * (1 - (1 - eps ^ 2) ^ 2 / 2 - sqrt(1 - (1 - eps ^ 2) ^
                                                    2)) / (1 - eps ^ 2) ^ 2
    AAu <-
      2 * AA0 * (1 - (1 - eps ^ 2) ^ 2 / 2 + sqrt(1 - (1 - eps ^ 2) ^
                                                    2)) / (1 - eps ^ 2) ^ 2
    return(c(AAl, AAu))
  }
  
  # HC
  obj_HC <- function(x, AA0, eps = grid.epsilon) {
    # Function for numerical search for the epsilon local grid for scaled half cauchy distributions
    return(integrate(function(t) {
      sqrt(dhalfcauchy(t, scale = x) * dhalfcauchy(t, scale = AA0))
    }, lower = 0, upper = Inf)$value - (1 - eps ^ 2))
  }
  
  # Lomax (LMX)
  obj_LMX <- function(x, AA0, eps = grid.epsilon) {
    # Function for numerical search for the epsilon local grid for scaled Lomax distributions
    return(integrate(function(t) {
      sqrt(dlomax(t, scale = x) * dlomax(t, scale = AA0))
    }, lower = 0, upper = Inf)$value - (1 - eps ^ 2))
  }
  
  parameters <- list()
  # tail_prob adjusting
  for (i in seq_along(distributions)) {
    if (!distributions[i] %in% c("HN", "HC", "EXP", "LMX"))
      stop(paste0(
        "The distribution ",
        distributions[i],
        " is not supported by this function."
      ))
    
    j <- (i - 1) * 2 + 1
    if (distributions[i] == "HN") {
      p_HN_l <- HN_A0_2_Al_Au(AA0 = AA0.HN)[1]
      p_HN_u <- HN_A0_2_Al_Au(AA0 = AA0.HN)[2]
      parameters <- c(parameters, p_HN_l, p_HN_u)
      names(parameters)[j:(j + 1)] <- c("p_HN_l", "p_HN_u")
    }
    if (distributions[i] == "HC") {
      p_HC_l <-
        uniroot(
          obj_HC,
          lower = 0.0001,
          upper = AA0.HC,
          tol = 1e-9,
          AA0 = AA0.HC,
          eps = grid.epsilon
        )$root
      p_HC_u <-
        uniroot(
          obj_HC,
          lower = AA0.HC,
          upper = 100,
          tol = 1e-9,
          AA0 = AA0.HC,
          eps = grid.epsilon
        )$root
      parameters <- c(parameters, p_HC_l, p_HC_u)
      names(parameters)[j:(j + 1)] <- c("p_HC_l", "p_HC_u")
    }
    if (distributions[i] == "EXP") {
      p_EXP_l <- EXP_A0_2_Al_Au(AA0 = AA0.EXP)[1]
      p_EXP_u <- EXP_A0_2_Al_Au(AA0 = AA0.EXP)[2]
      parameters <- c(parameters, p_EXP_l, p_EXP_u)
      names(parameters)[j:(j + 1)] <- c("p_EXP_l", "p_EXP_u")
    }
    if (distributions[i] == "LMX") {
      p_LMX_l <-
        uniroot(
          obj_LMX,
          lower = 0.0001,
          upper = AA0.LMX,
          tol = 1e-9,
          AA0 = AA0.LMX,
          eps = grid.epsilon
        )$root
      p_LMX_u <-
        uniroot(
          obj_LMX,
          lower = AA0.LMX,
          upper = 100,
          tol = 1e-9,
          AA0 = AA0.LMX,
          eps = grid.epsilon
        )$root
      parameters <- c(parameters, p_LMX_l, p_LMX_u)
      names(parameters)[j:(j + 1)] <- c("p_LMX_l", "p_LMX_u")
    }
  }
  
  return(parameters)
}