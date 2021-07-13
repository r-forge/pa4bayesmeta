\name{prior_accuracy}
\alias{prior_accuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sensitivity-based identification and classification of inaccurate heterogeneity priors
%%  ~~function to do ... ~~
}
\description{
This function implements the algorithm for sensitivity-based identification and classification of inaccurate heterogeneity priors 
proposed in Ott et al. (2021).
It can be applied to any proper heterogeneity prior which supports the assumption that small heterogeneity values are likely.
This function classifies the specified heterogeneity prior as either anticonservative (puts more than half of its probability mass on too small heterogeneity values) or conservative (puts more than half of its probability mass on too large heterogeneity values).
Assumes a Bayesian normal-normal hierarchical model.
}
\usage{
prior_accuracy(df, r.tau.prior, 
               MM = 10^6, mu.mean = 0, mu.sd = 4, 
               tail.prob = 0.5, grid.epsilon = 0.00354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "y" containing the effect estimates for the individual studies and a column "sigma" containing the standard errors of these estimates
}
\item{r.tau.prior}{
  function to generate random samples from the heterogeneity prior of interest
} 
\item{MM}{
  size of the random sample from the prior specified in \code{r.tau.prior}. Postitive integer. Defaults to 10^6.
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
First, it computes the effective median RLMC for the prior specified in \code{r.tau.prior} and the data set \code{df}.
Second, it applies the 50\%-RLMC-based adjustment to obtain a HN and a HC prior with the same median RLMC for the data set \code{df}.
Third, it computes epsilon-local sensitivity estimates for these two adjusted HN and HC priors and the data set \code{df}.
Finally, comparing the two sensitivity estimates leads to the decision:
If the sensitivity estimate for the HN prior is larger than the sensitivity estimate for the HC prior, then the specified prior is classfied as anticonservative. 
If the sensitivity estimate for the HC prior is larger than the sensitivity estimate for the HN prior, then the specified prior is classfied as conservative. 

The algorithm is described in Ott et al. (2021, Section 5), see the Supplementary Material of that manuscript (Section 4.6) for implementation details and an example.

If the heterogeneiy prior of interest is a HN or a HC prior, the algorithm described above can be simplified. In such a case, we recommend to use the \code{HN_accuracy} or the \code{HC_accuracy} function instead.
}
\value{
A list of the following 4 elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{param }{a list of the following two scale parameter values (non-negative real numbers): \cr
    \code{HN}: scale parameter of the adjusted HN prior \cr
    \code{HC}: scale parameter of the adjusted HC prior
  }
  \item{mrlmc_prior }{the effective median RLMC for the specified prior}
  \item{S_tau }{a list of the following two epsilon-local sensitivity estimates with respect to the heterogeneity parameter tau: \cr
    \code{`S(tau)^HN`}: sensitivity estimate for the HN prior \cr
    \code{`S(tau)^HC`}: sensitivity estimate for the HC prior
  }
  \item{decision }{a string specifiying the decision if the heterogeneity prior is anticonservative or conservative}
%  \item{p_HN }{scale parameter of half-normal prior}
%  \item{p_EXP }{scale parameter of exponential prior}
%  \item{p_HC }{scale parameter of half-Cauchy prior}
%  \item{p_LMX }{scale parameter of Lomax prior with shape parameter=1}
%% ...
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis.
In preparation. 2021.
%Manuscript submitted to \emph{Statistical Methods in Medical Research}. 


Ott, M., Hunanyan, S., Held, L., Roos, M. Supplementary Material: Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. In preparation. 2021.
%Submitted to \emph{Statistical Methods in Medical Research}. 

}
\seealso{
\code{\link{HN_accuracy}},
\code{\link{HC_accuracy}},
\code{\link{effective_rlmc}},
\code{\link{pri_par_adjust}}
%for the alternative static tail adjustment of priors
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load the steriod-resistant rejection (SRR) data analyzed in Friede et al. (2017)
data(srr)

# sensitivity-based accuracy estimation & classification for
# the uniform heterogeneity prior on [0, 4] and the SRR data
prior_accuracy(df = srr, 
               r.tau.prior = function(t) runif(t, min = 0, max = 4),
               mu.mean = 0, mu.sd = 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
