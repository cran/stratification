\name{strata.bh}
\alias{strata.bh}
\title{ Stratification of a Population Given a Set of Boundaries }
\description{
The function \code{strata.bh} stratifies a population given a set of boundaries. It calculates the stratum sample sizes and the anticipated coefficient of variation or relative root mean squared error.
}
\usage{
strata.bh(x, bh, n = NULL, CV = NULL, Ls = 3, certain = NULL,
          alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), takenone = 0, 
          bias.penalty = 1, takeall = 0, takeall.adjust = FALSE, 
          rh = rep(1, Ls), model = c("none", "loglinear", "linear",
          "random"), model.control = list())      
}
\arguments{
  \item{x}{ A vector containing the values of the stratification variable \eqn{X} for every unit in the population. }
  \item{bh}{ A vector of the \eqn{L-1} stratum boundaries where \eqn{L} is the total number of strata. Therefore, if \code{takenone=0} then \eqn{L}=\code{Ls}, and if \code{takenone=1} then \eqn{L}=\code{Ls}+1. }
  \item{n}{ A numeric: the target sample size. It has no default value. The argument \code{n} or the argument \code{CV} must be input. }
  \item{CV}{ A numeric: the target coefficient of variation or relative root mean squared error if \code{takenone}=1. It has no default value. The argument \code{CV} or the argument \code{n} must be input. }
  \item{Ls}{ A numeric: the number of sampled strata (excluding the take-none stratum). The default is 3. }
  \item{certain}{ A vector giving the position, in the vector \code{x}, of the units that must be included in the sample (see \code{\link{stratification-package}}). By default \code{certain} is \code{NULL}, which means that no units are a priori chosen to be in the sample. }
  \item{alloc}{ A list specifying the allocation scheme. The list must contain 3 numerics for the 3 exponents \code{q1}, \code{q2} and \code{q3} in the general allocation scheme (see \code{\link{stratification-package}}). The default is Neyman allocation (\code{q1}=\code{q3}=0.5 and \code{q2}=0)}
  \item{takenone}{ A numeric: the number of take-none strata (0 or 1). The default is 0, i.e. no take-none stratum is included.}
  \item{bias.penalty}{ A numeric between 0 and 1 giving the penalty for the bias in the anticipated mean squared error (MSE) of the survey estimator (see \code{\link{stratification-package}}). This argument is relevant only if \code{takenone}=1. The default is 1.}
  \item{takeall}{ A numeric: the number of take-all strata. The default is 0, i.e. no take-all stratum is included.}
  \item{takeall.adjust}{ A logical. If TRUE, when \eqn{n_h > N_h}{nh>Nh} for a take-some stratum, the \code{takeall} argument is increased by one and the stratification is carried out again. This is done as long as \eqn{n_h \leq N_h}{nh<=Nh} for every take-some stratum. The default is \code{FALSE}, i.e. no adjustment is made. Note: in other functions of the package \pkg{stratification}, this adjustment is not optional; it is made automatically (see \code{\link{stratification-package}}).}
  \item{rh}{ A vector giving the anticipated response rates in each of the \code{Ls} sampled strata. A single number can be given if the rates do not vary among strata. The default is 1 in each stratum. }
  \item{model}{ A character string identifying the model used to describe the discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can be \code{"none"} if one assumes \eqn{Y=X}, \code{"loglinear"} for the loglinear model with mortality, \code{"linear"} for the heteroscedastic linear model or \code{"random"} for the random replacement model (see \code{\link{stratification-package}} for a description of these models). The default is \code{"none"}.}
  \item{model.control}{ A list of model parameters (see \code{\link{stratification-package}}). The default values of the parameters correspond to the model \eqn{Y=X}. }
}
\value{
  \item{Nh }{A vector of length \eqn{L} containing the population sizes \eqn{N_h}{Nh}, i.e. the number of units in each stratum.}
  \item{nh }{A vector of length \eqn{L} containing the sample sizes \eqn{n_h}{nh}, i.e. the number of units to sample in each stratum. See \code{\link{stratification-package}} for information about the rounding used to get these values.}
  \item{n }{The total sample size (\code{sum(nh)}).}
  \item{nh.nonint }{A vector of length \eqn{L} containing the non-integer values of the sample sizes \eqn{n_h}{nh}.}
  \item{certain.info }{A vector giving statistics for the certain stratum (see \code{\link{stratification-package}}). It contains \code{Nc}, the number of units chosen a priori to be in the sample, and \code{meanc}, the anticipated mean of \eqn{Y} for these units. }
  \item{opti.nh }{The final value of the criteria to optimize : either the sum of the integer stratum sample sizes (if a target \code{CV} was given) or the RRMSE associated to integer stratum sample sizes (if a target \code{n} was given).}
  \item{opti.nhnonint }{The final value of the criteria to optimize : either the sum of the non-integer stratum sample sizes (if a target \code{CV} was given) or the RRMSE associated to non-integer stratum sample sizes (if a target \code{n} was given).}
  \item{meanh }{A vector of length \eqn{L} containing the anticipated means of \eqn{Y} in each stratum.}
  \item{varh }{A vector of length \eqn{L} containing the anticipated variances of \eqn{Y} in each stratum.}
  \item{mean }{A numeric: the anticipated global mean value of \eqn{Y}.}
  \item{RMSE }{A numeric: the root mean squared error (or standard error if \code{takenone}=0) of the anticipated global mean of \eqn{Y}. This is defined as the squared root of: (\code{bias.penalty} x bias of the mean)^2 + variance of the mean.}
  \item{RRMSE }{A numeric: the anticipated relative root mean squared error (or coefficient of variation if \code{takenone}=0) for the mean of \eqn{Y}, i.e. \code{RMSE} divided by \code{mean}. }
  \item{relativebias }{A numeric: the anticipated relative bias of the estimator, i.e. (\code{bias.penalty} x bias of the mean) divided by \code{mean}. If \code{takenone}=0, this numeric is zero. }
  \item{propbiasMSE }{A numeric: the proportion of the MSE attributable to the bias of the estimator, i.e. (\code{bias.penalty} x bias of the mean)^2 divided by the MSE of the \code{mean}. If \code{takenone}=0, this numeric is zero.}
  \item{stratumID }{A factor identifying, for each observation in the input \eqn{X}, the stratum to which it has been assigned. The value "certain" is given to units a priori chosen to be in the sample. }
  \item{takeall }{The number of take-all strata in the final solution. Note: It is possible that \eqn{n_h=N_h} for non take-all strata because the condition for an automatic addition of a take-all stratum is \eqn{n_h>N_h}. }
  \item{call }{The function call (object of class "call").}
  \item{date }{A character string that contains the system date and time when the function ended. }
  \item{args }{A list of all the arguments input to the function or used by default. }
}
\references{
Baillargeon, S. and  Rivest L.-P. (2011). The construction of stratified designs in R with the package stratification. \emph{Survey Methodology}, \bold{37}(1), 53-65. \url{ http://www.statcan.gc.ca/pub/12-001-x/2011001/article/11447-eng.pdf }
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{print.strata}}, \code{\link{plot.strata}}, \code{\link{strata.cumrootf}}, \code{\link{strata.geo}}, \code{\link{strata.LH}} }
\examples{
data(USbanks)
adjust<-strata.geo(x=USbanks,CV=0.01,Ls=4,alloc=c(0.35,0.35,0))
adjust
adjust$nh.nonint
noadjust<-strata.bh(x=USbanks,bh=adjust$bh,CV=0.01,Ls=4,
          alloc=c(0.35,0.35,0),takeall=0,takeall.adjust=FALSE)
noadjust
noadjust$nh.nonint
# without the adjustment for a take-all stratum, n is smaller than
# with the adjustment, but the target CV is not reached.
}
\keyword{ survey }
