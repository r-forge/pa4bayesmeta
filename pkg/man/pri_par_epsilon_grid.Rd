\name{pri_par_epsilon_grid}
\alias{pri_par_epsilon_grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grid computation for epsilon-local sensitivity
%%  ~~function to do ... ~~
}
\description{
Starting with a base prior from a one-parameter family, this function determines the parameters of two
priors from the same parametric family which have the given Hellinger distance to the base prior.
The supported one-parameter distributions are half-normal (HN), half-Cauchy (HC),
exponential (EXP) and Lomax (LMX) with shape parameter = 1.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pri_par_epsilon_grid(distributions = c("HN", "HC"),
                     AA0.HN = 1,  AA0.HC = 1, 
                     AA0.EXP = NULL, AA0.LMX = NULL, 
                     grid.epsilon = 0.00354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{distributions}{
 vector of strings. Specifies the parametric prior distributions to use. The options are "HN" (half-normal), "HC" (half-Cauchy), "EXP" (exponential) and "LMX" (Lomax with shape parameter = 1). Defaults to c("HN", "HC").
}
\item{AA0.HN}{
 scale parameter of the base half-normal prior
}
\item{AA0.EXP}{
  scale parameter of the base exponential prior
}
\item{AA0.HC}{
  scale parameter of the base half-Cauchy prior
}
\item{AA0.LMX}{
  scale parameter of the base Lomax prior with shape parameter=1
}
\item{grid.epsilon}{
  value for epsilon, the Hellinger distance between the base prior
  and the \cr
  "lower"/"upper" prior on the grid
}
}
\details{
This 1-dimensional epsilon grid is needed for computation of the
epsilon-local sensitivity in the functions \code{prior_accuracy}, \code{HN_accuracy} and \code{HC_accuracy}.
See Ott et al. (2021, Section 3.5) for a brief description of this methodology and the Supplementary Material of Ott et al. (2021, Section 4.5) and Roos et al. (2015) for more details.
The default value for \code{grid.epsilon} corresponds to the 
Hellinger distance between two normal distribution with unit variance and
a difference in means of 0.01.
See Roos et al. (2015, Section 2.2) and Ott et al. (2021b, Section 3.5) for calibration and interpretation of Hellinger distance values.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of maximum eight scale parameter values (including only those parameters for the distributions specified in \code{distributions}, in the same order as in \code{distributions}):
\item{p_HN_l }{parameter of the "lower" half-normal prior on the grid}
\item{p_HN_u }{parameter of the "upper" half-normal prior on the grid}
\item{p_HC_l }{parameter of the "lower" half-Cauchy prior on the grid}
\item{p_HC_u }{parameter of the "upper" half-Cauchy prior on the grid}
\item{p_EXP_l }{parameter of the "lower" exponential prior on the grid}
\item{p_EXP_u }{parameter of the "upper" exponential prior on the grid}
\item{p_LMX_l }{scale parameter of the "lower" Lomax prior (with shape parameter=1) on the grid}
\item{p_LMX_u }{scale parameter of the "upper" Lomax prior (with shape parameter=1) on the grid}
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. In preparation. 2021.
% Manuscript submitted to \emph{Statistical Methods in Medical Research}.

Ott, M., Hunanyan, S., Held, L., Roos, M. Supplementary Material: Sensitivity-based identification of inaccurate heterogeneity priors in Bayesian meta-analysis. In preparation. 2021.
%Submitted to \emph{Statistical Methods in Medical Research}.

Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/journals/bayesian-analysis/volume-10/issue-2/Sensitivity-Analysis-for-Bayesian-Hierarchical-Models/10.1214/14-BA909.full}

Ott, M., Plummer, M. and Roos M. (2021b). How vague is vague? How informative is informative?
Reference analysis for Bayesian meta-analysis. \emph{Statistics in Medicine} (accepted May 10, 2021). \url{https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9076}
%% ~put references to the literature/web site here ~
}

%\seealso{
%\code{\link{cal_h_dist}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# epsilon grid for the HN(0.5) and HC(1) base priors and 
# epsilon = 0.00354 (default)
pri_par_epsilon_grid(AA0.HN = 0.5, AA0.HC = 1)

# epsilon grid for the EXP(1) and LMX(2) (with shape parameter = 1) 
# base priors and epsilon = 0.00354 (default)
pri_par_epsilon_grid(distributions = c("EXP", "LMX"),
                     AA0.EXP = 1, AA0.LMX = 2)
                     
# load the steriod-resistant rejection (SRR) data analyzed in Friede et al. (2017)
data(srr)

# epsilon grid for 50\%-RLMC adjusted HN and HC base priors with 
# target RLMC = 0.5 and epsilon = 0.00354 (default) for SRR data
pri_par_epsilon_grid(AA0.HN = pri_par_adjust(df = srr)[[1]], 
                     AA0.HC = pri_par_adjust(df = srr)[[2]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

