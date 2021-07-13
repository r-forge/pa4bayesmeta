\name{sigma_ref}
\alias{sigma_ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reference standard deviation
%%  ~~function to do ... ~~
}
\description{
Computes the reference standard deviation of the
given meta-analysis data set.
Depending on the argument \code{type.sigma.ref}, either 
the geometric or a weighted harmonic mean is used.
%The reference standard deviation is defined as the geometric mean
%of the standard deviations of the individual studies
%(Sorbye & Rue 2014, equation (7)).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigma_ref(df, type.sigma.ref = "geometric")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "sigma" containing the standard errors of the estimates for the individual studies
}
  \item{type.sigma.ref}{
   either \code{"geometric"} or \code{"harmonic"}. Defaults to \code{"geometric"}.
   Specifies if the geometric mean or a weighted harmonic mean is used to compute
   the reference standard deviation.
   See details for more information.
}
}
\details{
  If \code{type.sigma.ref = "geometric"}, the reference standard deviation is given by the geometric mean
of the standard deviations of the individual studies
 (Sorbye & Rue 2014 (equation (7)). \cr
 If \code{type.sigma.ref = "harmonic"}, the reference standard deviation \eqn{\sigma_{ref}} is
 the square root of a weighted harmonic mean of the variances of the individual studies,
 as described in Hoaglin (2016, page 490).
 More precisely, we have
 \deqn{\sigma_{ref} = \sqrt{ (k-1) \sum w_i /((\sum w_i)^2 - \sum w_i^2)},}
 where \eqn{k} is the number of studies in the data frame \code{df} and the weights are
 \eqn{w_i=\sigma_i^{-2}}, \eqn{i =1, ... , k}, for the standard deviations \eqn{\sigma_i} (or standard errors)
 of the individual studies.
}
\value{
The reference standard deviation of the data set. Non-negative real number.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Sorbye, S., Rue, H. (2014). Scaling intrinsic Gaussian Markov random field priors in
spatial modelling. \emph{Spatial Statistics} \bold{8}, 39--51.
\url{https://doi.org/10.1016/j.spasta.2013.06.004}

Hoaglin, D. (2016). Misunderstandings about Q and "Cochran's Q test" in meta-analysis. 
\emph{Statistics in Medicine} \bold{35(4)}, 485--495.
\url{https://doi.org/10.1002/sim.6632}
%% ~put references to the literature/web site here ~
}

\examples{
# load the steriod-resistant rejection (SRR) data analyzed in Friede et al. (2017)
data(srr)

# geometric mean of with-study standard deviations in SRR data
sigma_ref(df = srr)
# harmonic mean of with-study standard deviations in SRR data
sigma_ref(df = srr, type.sigma.ref = "harmonic")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

