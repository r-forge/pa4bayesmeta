
HC_accuracy <- function(df,
                        scale.HC = 1,
                        mu.mean = 0,
                        mu.sd = 4,
                        tail.prob = 0.5,
                        grid.epsilon = 0.00354) {
  # sensitivity-based quantification of inaccurate heterogeneity prior assumptions in Bayesian meta-analysis
  # NNHM for the HC heterogeneity prior
  # input:
  # df: data frame in form compatible with bayesmeta
  # scale.HC: scale parameter of the HC heterogeneity prior
  # mu.mean: mean of the Normal prior for mu
  # mu.sd: sd of the Normal prior for mu (according to Roever(2018) Bayesmeta-Paper (unit-information prior for log-odds-ratios))
  # tail.prob=0.5: 50%-RLMC-adjustment of HN heterogeneity priors for tau with respect to the median RLMC induced by the HC prior
  # grid.epsilon=0.00354: grid.epsilon as suggested by Roos et al. (2015) corresponds to H(N(0,1),N(0.01,1))=0.00354
  # output:
  # list containing parameters of HN and HC, the effective median RLMC induced by the HC prior, S(tau)^HN and S(tau)^HC values
  # and the decision whether these priors are anti-conservative or conservative
  
  # checking of input values
  if(! scale.HC > 0)
    stop("The scale parameter scale.HC must be a positive real number.")
  
  if(! mu.sd > 0)
    stop("The standard deviation mu.sd must be a positive real number.")
  
  if(! (tail.prob >= 0 & tail.prob <= 1))
    stop("tail.prob must be in [0, 1].")
  
  if(! grid.epsilon > 0)
    stop("grid.epsilon must be a positive real number.")
  
 
  
  ####----  Adjusting of the HC heterogeneity prior for the RLMC induced by the HN prior ----####
  
  # compute the median RLMC induced by the HC(scale.HC) prior
  mrlmc_HC <- median(
    effective_rlmc(
      df = df,
      r.tau.prior = function(n)
        rhalfcauchy(n = n, scale = scale.HC),
      MM = 10 ^ 6,
      output = "sample"
    )
  )
  
  # adjust the scale parameter of the HC heterogeneity prior for the mrlmc_HN value
  tau_HN_dyn <- pri_par_adjust(
    df = df,
    distributions = "HN",
    rlmc = mrlmc_HC,
    tail.prob = tail.prob
  )[[1]]
  
  # save the scale parameters in a list
  scaling_param = list(HN = tau_HN_dyn, HC = scale.HC)
  
  
  ####---- Epsilon-local sensitivity for HN(scale_HN) and 50%-RLMC-adjusted HC(tau_HC_dyn) ----####
  
  ####---- Epsilon-local grid for A*|X| scaled distributions for tau  ----####
  
  #
  # grid.epsilon<-0.00354
  #
  
  # epsilon-local grid
  gr2p_HNHC_dyn <-
    pri_par_epsilon_grid(
      distributions = c("HN", "HC"),
      AA0.HN = scaling_param[["HN"]],
      AA0.HC = scaling_param[["HC"]],
      grid.epsilon = grid.epsilon
    )
  
  
  ####---- Epsilon-local sensitivity: light-tailed HN prior for tau ----####
  
  
  # HN_dyn:
  fit.bayesmeta.HN_dyn <- bayesmeta(
    y = df[, "y"],
    sigma = df[, "sigma"],
    mu.prior.mean = mu.mean,
    mu.prior.sd = mu.sd,
    tau.prior = function(t) {
      dhalfnormal(t, scale = scaling_param[["HN"]])
    }
  )
  
  fit.bayesmeta.HN_tau_l_dyn <- bayesmeta(
    y = df[, "y"],
    sigma = df[, "sigma"],
    mu.prior.mean = mu.mean,
    mu.prior.sd = mu.sd,
    tau.prior = function(t) {
      dhalfnormal(t, scale = gr2p_HNHC_dyn$p_HN_l)
    }
  )
  
  fit.bayesmeta.HN_tau_u_dyn <- bayesmeta(
    y = df[, "y"],
    sigma = df[, "sigma"],
    mu.prior.mean = mu.mean,
    mu.prior.sd = mu.sd,
    tau.prior = function(t) {
      dhalfnormal(t, scale = gr2p_HNHC_dyn$p_HN_u)
    }
  )
  
  integrand_HN_l_base_tau_dyn <-
    function(x) {
      sqrt(
        fit.bayesmeta.HN_dyn$dposterior(tau = x) * fit.bayesmeta.HN_tau_l_dyn$dposterior(tau =
                                                                                           x)
      )
    }
  sens_HN_l_base_tau_dyn <-
    sqrt(1 - integrate(
      integrand_HN_l_base_tau_dyn,
      lower = 0,
      upper = Inf
    )$value) / grid.epsilon
  integrand_HN_u_base_tau_dyn <-
    function(x) {
      sqrt(
        fit.bayesmeta.HN_dyn$dposterior(tau = x) * fit.bayesmeta.HN_tau_u_dyn$dposterior(tau =
                                                                                           x)
      )
    }
  sens_HN_u_base_tau_dyn <-
    sqrt(1 - integrate(
      integrand_HN_u_base_tau_dyn,
      lower = 0,
      upper = Inf
    )$value) / grid.epsilon
  worst_sens_HN_tau_dyn <-
    max(sens_HN_l_base_tau_dyn, sens_HN_u_base_tau_dyn)
  
  
  ####---- Epsilon-local sensitivity: heavy-tailed HC prior for tau ----####
  
  # HC_dyn:
  fit.bayesmeta.HC_dyn <- bayesmeta(
    y = df[, "y"],
    sigma = df[, "sigma"],
    mu.prior.mean = mu.mean,
    mu.prior.sd = mu.sd,
    tau.prior = function(t) {
      dhalfcauchy(t, scale = scaling_param[["HC"]])
    }
  )
  
  fit.bayesmeta.HC_tau_l_dyn <- bayesmeta(
    y = df[, "y"],
    sigma = df[, "sigma"],
    mu.prior.mean = mu.mean,
    mu.prior.sd = mu.sd,
    tau.prior = function(t) {
      dhalfcauchy(t, scale = gr2p_HNHC_dyn$p_HC_l)
    }
  )
  
  fit.bayesmeta.HC_tau_u_dyn <- bayesmeta(
    y = df[, "y"],
    sigma = df[, "sigma"],
    mu.prior.mean = mu.mean,
    mu.prior.sd = mu.sd,
    tau.prior = function(t) {
      dhalfcauchy(t, scale = gr2p_HNHC_dyn$p_HC_u)
    }
  )
  
  integrand_HC_l_base_tau_dyn <-
    function(x) {
      sqrt(
        fit.bayesmeta.HC_dyn$dposterior(tau = x) * fit.bayesmeta.HC_tau_l_dyn$dposterior(tau =
                                                                                           x)
      )
    }
  sens_HC_l_base_tau_dyn <-
    sqrt(1 - integrate(
      integrand_HC_l_base_tau_dyn,
      lower = 0,
      upper = Inf
    )$value) / grid.epsilon
  integrand_HC_u_base_tau_dyn <-
    function(x) {
      sqrt(
        fit.bayesmeta.HC_dyn$dposterior(tau = x) * fit.bayesmeta.HC_tau_u_dyn$dposterior(tau =
                                                                                           x)
      )
    }
  sens_HC_u_base_tau_dyn <-
    sqrt(1 - integrate(
      integrand_HC_u_base_tau_dyn,
      lower = 0,
      upper = Inf
    )$value) / grid.epsilon
  worst_sens_HC_tau_dyn <-
    max(sens_HC_l_base_tau_dyn, sens_HC_u_base_tau_dyn)
  
  
  
  ####---- put everything together ----####
  
  # save the sensitivity values in list and add the names of the elements
  sensitivity <- list(worst_sens_HN_tau_dyn,
                      worst_sens_HC_tau_dyn)
  names(sensitivity) <- c("S(tau)^HN", "S(tau)^HC")
  
  
  ####---- generate desision ----####
  if (sensitivity[["S(tau)^HN"]] > sensitivity[["S(tau)^HC"]]) {
    decision <- "S(tau)^HN > S(tau)^HC: HN and HC are anti-conservative"
  } else {
    decision <- "S(tau)^HN < S(tau)^HC: HN and HC are conservative"
  }
  
  return(
    list(
      param = scaling_param,
      mrlmc_HC = mrlmc_HC,
      S_tau = sensitivity,
      decision = decision
    )
  )
}