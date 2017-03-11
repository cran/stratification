\name{var.strata}
\alias{var.strata}
\alias{print.var.strata}
\title{ Anticipated Variances and RRMSE from a Stratified Design for a Survey Variable Y }
\description{
\code{var.strata} calculates the anticipated means, variances and relative root mean squared error (RRMSE) obtained when applying a stratified design to a survey variable \eqn{Y}. The variable \eqn{Y} can be input or it can be defined from \eqn{X} by a specified loglinear with mortality, heteroscedastic linear or random replacement model.\cr

\code{print.var.strata} prints a "var.strata" object, presenting the stratification information into a table. 
}
\usage{
var.strata(strata, y = NULL, rh = strata$args$rh, rh.postcorr = 
           FALSE, model = c("none", "loglinear", "linear", "random"), 
           model.control = list())

\method{print}{var.strata}(x, \dots)
}
\arguments{
  \item{strata}{ An object of class "strata", which represents a stratified design. }
  \item{y}{ A vector containing the values of the survey variable \eqn{Y} for every unit of the population, respecting the order of the units in the \code{x}-vector used to create \code{strata}. The default is that \eqn{Y} is not given. }
  \item{rh}{ A vector giving the anticipated response rates in each of the \code{Ls} sampled strata. A single number can be given if the rates do not vary among strata. The default is to use the rates given in the \code{strata.bh} object. }
  \item{rh.postcorr}{ A logical. If TRUE, a posterior correction for non-response is applied. This correction takes into account the non-response in the \code{strata.bh} object. It is only available when the stratified design \code{strata} had a target CV. The default is FALSE, i.e. no posterior correction is made (see \bold{Details}). }
  \item{model}{ A character string identifying the model used to describe the discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can be \code{"none"} if one assumes \eqn{Y=X}, \code{"loglinear"} for the loglinear model with mortality, \code{"linear"} for the heteroscedastic linear model or \code{"random"} for the random replacement model (see \code{\link{stratification-package}} for a description of these models). The default is \code{"none"}.}
  \item{model.control}{ A list of model parameters (see \code{\link{stratification-package}}). The default values of the parameters correspond to the model \eqn{Y=X}. }
  \item{x}{ An object of class "var.strata" to print. }
  \item{\dots}{ Additional arguments affecting the print produced. }
}
\details{
POSTERIOR CORRECTION FOR NON-RESPONSE (with a target CV only\cr

The optional posterior correction for non-response is done as follows. For each take-some stratum, \eqn{n_h}{nh} is increased if the input \code{rh} is lower than the anticipated response rate in the \code{strata.bh} object, and \eqn{n_h}{nh} is decreased if the input \code{rh} is higher than the anticipated response rate given when creating the \code{strata.bh} object. The modification of \eqn{n_h}{nh} is done by multiplying it by \code{strata$args$rh/rh}.\cr

The weakness of this posterior correction is that it cannot take into account non-response in a take-all stratum. In that stratum, \eqn{n_h}{nh} cannot be increased since it is equal to \eqn{N_h}{Nh}. To correctly account for non-response in a take-all stratum, the boundary of the stratum has to be lowered. This is what the generalized Lavallee-Hidiroglou method does (\code{\link{strata.LH}}).
}
\value{
  \item{nh }{A vector of length \code{L} containing the integer sample sizes \eqn{n_h}{nh}, i.e. the number of units to sample in each stratum. This vector can be different than \code{strata$nh} if \code{rh.postcorr=TRUE}.}
  \item{n }{The total sample size (\code{sum(nh)}). This number can be different than \code{strata$n} if \code{rh.postcorr=TRUE}.}
  \item{nhnonint }{A vector of length \eqn{L} containing the non-integer values of the sample sizes. This vector can be different than \code{strata$nhnonint} if \code{rh.postcorr=TRUE}.}
  \item{certain.info }{A vector giving statistics for the certainty stratum (see \code{\link{stratification-package}}). It contains \code{Nc}, the number of units chosen a priori to be in the sample, and \code{meanc}, the anticipated mean of \eqn{Y} for these units. }
  \item{meanh }{A vector of length \eqn{L} containing the anticipated means of \eqn{Y} in each stratum.}
  \item{varh }{A vector of length \eqn{L} containing the anticipated variances of \eqn{Y} in each stratum.}
  \item{mean }{A numeric: the anticipated global mean value of \eqn{Y}.}
  \item{RMSE }{A numeric: the root mean squared error (or standard error if \code{strata$args$takenone}=0) of the anticipated global mean of \eqn{Y}. This is defined as the squared root of: (\code{bias.penalty} x bias of the mean)^2 + variance of the mean.}
  \item{RRMSE }{A numeric: the anticipated relative root mean squared error (or coefficient of variation if  \code{strata$args$takenone}=0) for the mean of \eqn{Y}, i.e. \code{RMSE} divided by \code{mean}. }
  \item{relativebias }{A numeric: the anticipated relative bias of the estimator, i.e. (\code{bias.penalty} x bias of the mean) divided by \code{mean}. If  \code{strata$args$takenone}=0, this numeric is zero. }
  \item{propbiasMSE }{A numeric: the proportion of the MSE attributable to the bias of the estimator, i.e. (\code{bias.penalty} x bias of the mean)^2 divided by the MSE of the \code{mean}. If  \code{strata$args$takenone}=0, this numeric is zero.}
  \item{call }{The function call (object of class "call").}
  \item{date }{A character string that contains the system date and time when the function ended. }
  \item{args }{A list of all the arguments input to the function or used by default. }
}
\references{
Baillargeon, S. and  Rivest L.-P. (2011). The construction of stratified designs in R with the package stratification. \emph{Survey Methodology}, \bold{37}(1), 53-65.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{strata.bh}}, \code{\link{strata.cumrootf}}, \code{\link{strata.geo}}, \code{\link{strata.LH}} }
\examples{
nomodel <- strata.LH(x=Sweden$REV84, CV=0.05, Ls=3, alloc=c(0.5,0,0.5),
          takeall=1, model="none")
# We can give a vector of the Y values for every unit in the population
var.strata(nomodel, y=Sweden$RMT85) 
# Or specify a model between X and Y
var.strata(nomodel, model="loglinear", model.control=list(beta=1.058355,
           sig2=0.06593083, ph=1))
# Compared to taking into account the model in the optimization
model <- strata.LH(x=Sweden$REV84, CV=0.05, Ls=3, alloc=c(0.5,0,0.5),
         takeall=1, model="loglinear", model.control=list(beta=1.058355,
		     sig2=0.06593083, ph=1))
var.strata(model, y=Sweden$RMT85)

### Examples of posterior correction for non-response
LH <- strata.LH(x=MRTS, CV=0.01, Ls=4, alloc=c(0.5,0,0.5), takeall=1)
LH
# Without non-response in the take-all strata
var.strata(LH, rh.postcorr=TRUE, rh=c(0.85,0.9,0.9,1))
strata.LH(x=MRTS, CV=0.01, Ls=4, alloc=c(0.5,0,0.5), takeall=1, rh=c(0.85,0.9,0.9,1))
# With non-response in the take-all strata
var.strata(LH, rh.postcorr=TRUE, rh=0.9)
strata.LH(x=MRTS, CV=0.01, Ls=4, alloc=c(0.5,0,0.5), takeall=1, rh=0.9)
}
\keyword{ survey }
