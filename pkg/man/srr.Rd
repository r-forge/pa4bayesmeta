\name{srr}
\alias{srr}
\docType{data}
\title{
Steroid-resistant rejection data
%%   ~~ data name/kind ... ~~
}
\description{
Meta-analysis data set including 3 controlled studies comparing steroid-resistant rejections
under Interleukin-2
receptor antibodies (IL-2RA, Basiliximab or Daclizumab) treatment versus other immunosuppressant drugs after liver transplantation in children. The effect estimates of interest are (log) odds ratios.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(srr)}
\format{
  A data frame with 3 observations on the following 7 variables:
  \describe{
    \item{\code{study}}{character string, label of the study}
    %TODO: specify the direction of the OR
    \item{\code{n_events_t}}{integer, number of steroid-resistant rejections in the group treated with Interleukin-2 receptor antibodies}
    \item{\code{n_t}}{integer, number of children treated with Interleukin-2 receptor antibodies}
    \item{\code{n_events_c}}{integer, number of steroid-resistant rejections in the group treated with other immunosuppressant drugs}
    \item{\code{n_c}}{integer, number of children treated with other immunosuppressant drugs (control group)}
    \item{\code{y}}{numeric, log odds ratio estimate for the study}
    \item{\code{sigma}}{numeric, standard error of the log odds ratio estimate for the study}
  }
}
\details{
This data set was originally analyzed by Crins et al. (2014) and reanalyzed by Friede et al. (2017).
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Friede, T., Roever, C., Wandel, S., Neuenschwander, B. (2017). Meta-analysis of few
small studies in orphan diseases. \emph{Research Synthesis Methods} \bold{8(1)}, 79--91.
\url{https://doi.org/10.1002/jrsm.1217}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Crins N. D., Roever C., Goralczyk A. D., Friede T. (2014). Interleukin-2 receptor antagonists 
for pediatric liver transplant recipients: a systematic review and meta-analysis of controlled studies. \emph{Pediatric Transplantation} \bold{18}: 839--850.
\url{https://onlinelibrary.wiley.com/doi/abs/10.1111/petr.12362}
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(srr)
str(srr)
# forest plot
forest(x = srr$y, sei = srr$sigma, xlab = "log odds ratio")
}
\keyword{datasets}
