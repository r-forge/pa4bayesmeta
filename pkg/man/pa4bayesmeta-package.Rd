\name{pa4bayesmeta-package}
\alias{pa4bayesmeta-package}
\docType{package}
\title{
\packageTitle{pa4bayesmeta}
}
\description{
This package provides functions to identify and classify inaccurate heterogeneity priors (i.e. priors for the between-study standard deviation) in the Bayesian normal-normal hierarchical model used for Bayesian meta-analysis. It implements the methodology proposed in Ott et al. (2021).

A heterogeneity prior is called inaccurate if it does not assign equal probability mass to both sides of the true between-study standard deviation. Whereas an anticonservative heterogeneity prior puts more than half of its probability mass on heterogeneity values that are smaller than the true value, a conservative heterogeneity prior puts more than half of its probability mass on heterogeneity values that are larger than the true value.

The main function \code{\link{prior_accuracy}}() provides sensitivity estimates for identification of inaccurate heterogeneity priors and classifies the specified heterogeneity prior as either anticonservative or conservative.
For the special cases of half-normal (HN) and half-Cauchy (HC) heterogeneity priors, the functions \code{\link{HN_accuracy}}() and \code{\link{HN_accuracy}}() are recommended, which apply a simplified algorithm.

The function \code{\link{effective_rlmc}}() computes the relative latent model complexity (RLMC) associated with a heterogeneity prior and a data set.
The function \code{\link{pri_par_adjust}}() implements the novel 50\%-RLMC-based prior adjustment introduced in Ott et al. (2021). This function determines the scale parameter of certain one-parameter distributions (HN, HC, EXP (exponential), LMX (Lomax) with shape parameter = 1) for a given target RLMC, a tail probability and a data set.
In order to unify notation, the heterogeneity priors are defined
as scaled distributions \eqn{\tau \sim A_0 |X|}{tau ~ A_0 |X|},
where \eqn{A_0} is a scale parameter and
\eqn{X} is the standard form of the distribution.

The package also contains two medical meta-analysis data sets including a small (2 or 3) number of studies, which are used in Ott et al. (2021):
the steroid-resistant rejection (SRR) and the kidney disease (KD) data sets.
These data sets can be loaded by specifying \code{data(srr)} and \code{data(kd)}.
}
\details{
Package: pa4bayesmeta

Type: Package

Title: Prior accuracy for Bayesian meta-analysis

Version: 0.1-4

Date: 2021-08-01

Author: Manuela Ott [aut, cre], Malgorzata Roos [aut] 

Maintainer: Manuela Ott <manuela.c.ott@gmail.com>

Depends: bayesmeta

License: GPL (>=2)
%The DESCRIPTION file:
%\packageDESCRIPTION{pCalibrate}
%\packageIndices{pCalibrate}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Manuela Ott, Malgorzata Roos
%\packageAuthor{pCalibrate}
Maintainer: \packageMaintainer{pa4bayesmeta}
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. 
Submitted to Statistical Methods in Medical Research. 2021.
%Manuscript submitted to \emph{Statistical Methods in Medical Research}. 

Ott, M., Hunanyan, S., Held, L., Roos, M. Supplementary Material: Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. Submitted to Statistical Methods in Medical Research. 2021.

%Spiegelhalter, D., Abrams, K., Myles, J. (2004). Bayesian Approaches to Clinical Trials and Health-Care Evaluation. John Wiley & Sons, Ltd.
}
\keyword{package}
%\keyword{htest}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# load the steriod-resistant rejection (SRR) data analyzed in Friede et al. (2017)
data(srr)

# sensitivity-based accuracy estimation & classification for
# the uniform heterogeneity prior on [0, 4] and the SRR data
prior_accuracy(df = srr, 
               r.tau.prior = function(t) runif(t, min = 0, max = 4))
               
# sensitivity-based accuracy estimation & classification for
# the HN(1) prior and the SRR data
HN_accuracy(df = srr, scale.HN = 1,
            mu.mean = 0, mu.sd = 4)
            
# summary statistics for MC sample of RLMC values
# for the HN(0.5) prior and the SRR data
effective_rlmc(df = srr, r.tau.prior = function(n) rhalfnormal(n = n, scale = 0.5), 
               output = "summary")
               
# 50\%-RLMC-based adjustment of HN and HC priors used in Ott et al. (2021)
# with target RLMC 0.5
pri_par_adjust(df = srr, rlmc = 0.5)
}
