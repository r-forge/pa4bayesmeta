\name{pri_par_adjust}
\alias{pri_par_adjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heterogeneity prior adjustment based on the relative latent model complexity (RLMC)
%%  ~~function to do ... ~~
}
\description{
Adjusts the scale parameter of 
the specified one-parameter priors
for the between-study standard deviation
such that the probability mass above the reference threshold
equals the given tail probability.
The reference threshold used depends on the target RLMC and the with-study
standard errors in the data set.
The supported priors are the half-normal (HN), the half-Cauchy (HC),
the exponential (EXP) and the Lomax (LMX) distributions.
The shape parameter of the LMX distribution is fixed at 1.
}
\usage{
pri_par_adjust(df, rlmc = 0.5, tail.prob = 0.5, 
               distributions = c("HN", "HC"), 
               type.sigma.ref = "geometric")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "sigma" containing the standard errors of the estimates for the individual studies
}
\item{rlmc}{
  target relative latent model complexity. Real number in [0,1]. Defaults to 0.5.
}
\item{tail.prob}{
  probability mass of the prior above the reference threshold (which depends on \code{rlmc} and \code{df$sigma}) Real number in [0,1]. Defaults to 0.5, which corresponds to aliging the median of the prior.
}
\item{distributions}{
 vector of strings. Specifies the parametric prior distributions to use. The options are "HN" (half-normal), "HC" (half-Cauchy), "EXP" (exponential) and "LMX" (Lomax with shape parameter = 1). Defaults to c("HN", "HC").
}
\item{type.sigma.ref}{
   either \code{"geometric"} or \code{"harmonic"}. Defaults to \code{"geometric"}.
   Specifies if the geometric mean or a weighted harmonic mean is used to compute
   the reference standard deviation.
   See details for more information.
}
}
\details{
This heterogeneity prior adjustment applies to Bayesian meta-analysis
expressed by a normal-normal hierarchical model.
The reference threshold U is given by
\deqn{U = \sigma_{ref} \sqrt{rlmc/(1-rlmc)},}
where \eqn{\sigma_{ref}} is the reference standard deviation of the data set,
i.e. by default the geometric mean of \code{df$sigma}.
Then, a prior tail-adjustment is applied for this reference threshold
and the specified tail probability:
The free parameter of the prior is determined such that
\deqn{P[\tau > U] = \alpha,}
where \eqn{\alpha =} \code{tail.prob}.

Ott et al. (2021) suggest to use \code{tail.prob = 0.5} as default,
so that the median of the prior will be aligned with the reference threshold.
%%  ~~ If necessary, more details than the description above ~~

If \code{type.sigma.ref = "geometric"}, the reference standard deviation is given by the geometric mean
of the standard deviations of the individual studies
 (Sorbye & Rue 2014 (equation (7)). \cr
 If \code{type.sigma.ref = "harmonic"}, the reference standard deviation \eqn{\sigma_{ref}} is
 the square root of a weighted harmonic mean of the variances of the individual studies,
 as described in Hoaglin (2016, page 490).
 See \code{\link{sigma_ref}} for the formula.
}
\value{
A list of maximum four scale parameter values (including only those parameters for the distributions specified in \code{distributions}, in the same order as in \code{distributions}):
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{p_HN }{parameter of half-normal prior}
 \item{p_HC }{parameter of half-Cauchy prior}
 \item{p_EXP }{parameter of exponential prior}
 \item{p_LMX }{scale parameter for Lomax prior with shape parameter=1}
%  \item{p_HN }{scale parameter of half-normal prior}
%  \item{p_EXP }{scale parameter of exponential prior}
%  \item{p_HC }{scale parameter of half-Cauchy prior}
%  \item{p_LMX }{scale parameter of Lomax prior with shape parameter=1}
%% ...
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. In preparation. 2021.
%Manuscript submitted to \emph{Statistical Methods in Medical Research}.

Sorbye, S., Rue, H. (2014). Scaling intrinsic Gaussian Markov random field priors in
spatial modelling. \emph{Spatial Statistics} \bold{8}, 39--51.
\url{https://doi.org/10.1016/j.spasta.2013.06.004}

Hoaglin, D. (2016). Misunderstandings about Q and "Cochran's Q test" in meta-analysis. 
\emph{Statistics in Medicine} \bold{35(4)}, 485--495.
\url{https://doi.org/10.1002/sim.6632}
}

\seealso{
\code{\link{sigma_ref}}
%for the alternative static tail adjustment of priors
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load the steriod-resistant rejection (SRR) data analyzed in Friede et al. (2017)
data(srr)
 
# for the SRR data, compute the following:                 
# 50\%-RLMC-based adjustment of HN and HC priors used in Ott et al. (2021)
# with target RLMC 0.5
pri_par_adjust(df = srr, rlmc = 0.5)
# 50\%-RLMC-based adjustment of EXP and LMX priors used in Ott et al. (2021)
# with target RLMC 0.8
pri_par_adjust(df = srr, distributions = c("EXP", "LMX"),
               rlmc = 0.8)
# 50\%-RLMC-based adjustment of HN and HC priors with target RLMC 0.2
# using the harmonic mean to determine the reference threshold
pri_par_adjust(df = srr, rlmc = 0.2,
              type.sigma.ref = "harmonic")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
