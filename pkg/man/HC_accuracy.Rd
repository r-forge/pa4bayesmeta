\name{HC_accuracy}
\alias{HC_accuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sensitivity-based identification and classification of inaccurate half-Cauchy heterogeneity priors
%%  ~~function to do ... ~~
}
\description{
This function implements the algorithm for sensitivity-based identification and classification of inaccurate half-Cauchy heterogeneity priors 
proposed in Ott et al. (2021).
It classifies the specified half-Cauchy heterogeneity prior as either anticonservative (puts more than half of its probability mass on too small heterogeneity values) or conservative (puts more than half of its probability mass on too large heterogeneity values).
Assumes a Bayesian normal-normal hierarchical model.
}
\usage{
HC_accuracy(df, scale.HC = 1, mu.mean = 0, mu.sd = 4, 
            tail.prob = 0.5, grid.epsilon = 0.00354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "y" containing the effect estimates for the individual studies and a column "sigma" containing the standard errors of these estimates
}
\item{scale.HC}{
  scale parameter of the half-Cauchy prior of interest
} 
\item{mu.mean}{
  mean of the normal prior for the effect mu. Defaults to 0.
} 
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu. Defaults to 4.
} 
\item{tail.prob}{
  probability mass of the prior above the reference threshold. Defaults to 0.5, which corresponds to the 50\%-RLMC-based adjustment of priors.
} 
\item{grid.epsilon}{
  value for epsilon, a tuning parameter for epsilon-local sensitivity computation (see \code{\link{pri_par_epsilon_grid}} for details). Defaults to 0.00354.
}
}
\details{
This algorithm proceeds as follows:
First, it computes the effective median RLMC for the HC prior with scale parameter \code{scale.HC} and the data set \code{df}.
Second, it applies the 50\%-RLMC-based adjustment to obtain a HN prior with the same median RLMC for the data set \code{df}.
Third, it computes epsilon-local sensitivity estimates for these two HN and HC priors and the data set \code{df}.
Finally, comparing the two sensitivity estimates leads to the decision:
If the sensitivity estimate for the HN prior is larger than the sensitivity estimate for the HC prior, then both the HN and HC prior are classfied as anticonservative. 
If the sensitivity estimate for the HC prior is larger than the sensitivity estimate for the HN prior, then both the HN and HC prior are classfied as conservative. 

This is a special case of the algorithm described in Ott et al. (2021, Section 5), see also the Supplementary Material of that manuscript (Section 4.6) for more information.

Using HC heterogeneity priors as "weakly informative" priors has been suggested in Gelman (2006). The standard half-Cauchy heterogeneity prior (\code{scale.HC = 1}) has for example been used in Bodnar et al. (2017).
}
\value{
A list of the following 4 elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{param }{a list of the following two scale parameter values (non-negative real numbers): \cr
    \code{HN}: scale parameter of the adjusted HN prior \cr
    \code{HC}: scale parameter of the specified HC prior
    }
  \item{mrlmc_HC }{the effective median RLMC for the specified HC prior}
  \item{S_tau }{a list of the following two epsilon-local sensitivity estimates with respect to the heterogeneity parameter tau: \cr
    \code{`S(tau)^HN`}: sensitivity estimate for the HN prior \cr
    \code{`S(tau)^HC`}: sensitivity estimate for the HC prior
    }
  \item{decision }{a string specifiying the decision if the HN and HC priors are anticonservative or conservative}
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis.
In preparation. 2021.
%Manuscript submitted to \emph{Statistical Methods in Medical Research}. 

Ott, M., Hunanyan, S., Held, L., Roos, M. Supplementary Material: Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. In preparation. 2021.
%Submitted to \emph{Statistical Methods in Medical Research}. 

Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models (Comment on article by Browne and Draper). \emph{Bayesian Analysis} \bold{1(3)}, 515--534.
\url{https://doi.org/10.1214/06-BA117A}

Bodnar, O., Link, A., Arendacka, B., Possolo, A., Elster C. (2017). Bayesian estimation in random effects meta-analysis using a non-informative prior. \emph{Statistics in Medicine} \bold{36(2)}, 378--399. 
\url{https://doi.org/10.1002/sim.7156}
}
\seealso{
\code{\link{HN_accuracy}},
\code{\link{prior_accuracy}},
\code{\link{effective_rlmc}},
\code{\link{pri_par_adjust}}
%for the alternative static tail adjustment of priors
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load the steriod-resistant rejection (SRR) data analyzed in Friede et al. (2017)
data(srr)

# sensitivity-based accuracy estimation & classification for
# the HC(1) prior and the SRR data
HC_accuracy(df = srr, scale.HC = 1,
            mu.mean = 0, mu.sd = 4)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
}