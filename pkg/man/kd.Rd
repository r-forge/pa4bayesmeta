\name{kd}
\alias{kd}
\docType{data}
\title{
Kidney disease data
%%   ~~ data name/kind ... ~~
}
\description{
Meta-analysis data set including 2 controlled studies comparing Belatacept versus Ciclosporin A
for prophylaxis of graft rejection in adults receiving a renal transplant. The effect estimates of interest are (log) hazard ratios.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(kd)}
\format{
  A data frame with 2 observations on the following 3 variables:
  \describe{
    \item{\code{study}}{character string, label of the study}
    \item{\code{y}}{numeric, log hazard ratio estimate for the study}
    \item{\code{sigma}}{numeric, standard error of the log hazard ratio estimate for of the study}
  }
}
\details{
This data set was originally analyzed by the Institute for Quality and Efficiency in Health Care (2015) and reanalyzed by Bender et al. (2018).

Both studies included in this data set provide evidence for a benefit of Belatacept compared to Ciclosporin A.
In the original analysis, no pooling of the study results was done due to large heterogeneity.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Bender, R., Friede, T., Koch, A., Kuss, O., Schlattmann, P., Schwarzer, G., Skipka,
G. (2018). Methods for evidence synthesis in the case of very few studies. \emph{Research
Synthesis Methods} \bold{9(3)}, 382--392.
\url{https://doi.org/10.1002/jrsm.1297}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Institute for Quality and Efficiency in Health Care (2015). Belatacept -- Nutzenbewertung gemaess Paragraph 35a SGB V, Auftrag A15-25, Version 1.0 vom 13.10.2015.
\url{https://www.iqwig.de/projekte/a15-25.html}
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(kd)
str(kd)
# forest plot
forest(x = kd$y, sei = kd$sigma, xlab = "log hazard ratio")
}
\keyword{datasets}
