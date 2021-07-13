\name{HN_accuracy}
\alias{HN_accuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sensitivity-based identification and classification of inaccurate half-normal heterogeneity priors
%%  ~~function to do ... ~~
}
\description{
This function implements the algorithm for sensitivity-based identification and classification of inaccurate half-normal heterogeneity priors 
proposed in Ott et al. (2021).
It classifies the specified half-normal heterogeneity prior as either anticonservative (puts more than half of its probability mass on too small heterogeneity values) or conservative (puts more than half of its probability mass on too large heterogeneity values).
Assumes a Bayesian normal-normal hierarchical model.
}
\usage{
HN_accuracy(df, scale.HN = 0.5, mu.mean = 0, mu.sd = 4, 
            tail.prob = 0.5, grid.epsilon = 0.00354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "y" containing the effect estimates for the individual studies and a column "sigma" containing the standard errors of these estimates
}
\item{scale.HN}{
  scale parameter of the half-normal prior of interest. Defaults to 0.5.
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
First, it computes the effective median RLMC for the half-normal prior with scale parameter \code{scale.HN} and the data set \code{df}.
Second, it applies the 50\%-RLMC-based adjustment to obtain a HC prior with the same median RLMC for the data set \code{df}.
Third, it computes epsilon-local sensitivity estimates for these two HN and HC priors and the data set \code{df}.
Finally, comparing the two sensitivity estimates leads to the decision:
If the sensitivity estimate for the HN prior is larger than the sensitivity estimate for the HC prior, then both the HN and HC prior are classfied as anticonservative. 
If the sensitivity estimate for the HC prior is larger than the sensitivity estimate for the HN prior, then both the HN and HC prior are classfied as conservative. 

This is a special case of the algorithm described in Ott et al. (2021, Section 5), see also the Supplementary Material of that manuscript (Section 4.6) for more information.

The default value 0.5 of \code{scale.HN}, the scale parameter of the HN heterogeneity prior,
is recommended in Bender et al. (2018) for the meta-analysis of a small number of small studies.
}
\value{
A list of the following 4 elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{param }{a list of the following two scale parameter values (non-negative real numbers): \cr
    \code{HN}: scale parameter of the specified HN prior \cr
    \code{HC}: scale parameter of the adjusted HC prior
    }
  \item{mrlmc_HN }{the effective median RLMC for the specified HN prior. Real number in [0, 1].}
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

Bender, R., Friede, T., Koch, A., Kuss, O., Schlattmann, P., Schwarzer, G., Skipka,
G. (2018). Methods for evidence synthesis in the case of very few studies. \emph{Research
Synthesis Methods} \bold{9(3))}, 382--392.
\url{https://doi.org/10.1002/jrsm.1297}
}
\seealso{
\code{\link{HC_accuracy}},
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
# the HN(0.5) prior and the SRR data
HN_accuracy(df = srr, scale.HN = 0.5,
            mu.mean = 0, mu.sd = 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
