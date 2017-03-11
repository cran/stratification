\name{strata.rule}
\alias{strata.cumrootf}
\alias{strata.geo}
\title{ Non-Iterative Methods of Strata Construction }
\description{
These functions first determine boundaries to stratify a population. Then, in a second independent step, the sample sizes are calculated given a CV or the CV is computed given the total sample size. The function \code{strata.cumrootf} uses the cumulative root frequency method by Dalenius and Hodges (1959) and \code{strata.geo} uses the geometric method by Gunning and Horgan (2004). A model can be specified for the relationship between the stratification variable \eqn{X} and the survey variable \eqn{Y}, but this model has no impact on the first step of boundary determination. It only influences the calculation of the n or of the CV by the use of anticipated means and variances of \eqn{Y} instead of the empirical means and variances of \eqn{X}.
}
\usage{
strata.cumrootf(x, n = NULL, CV = NULL, Ls = 3, certain = NULL, 
       alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), rh = rep(1, Ls), 
       model = c("none", "loglinear", "linear", "random"),
       model.control = list(), nclass = NULL)

strata.geo(x, n = NULL, CV = NULL, Ls = 3, certain=NULL,
       alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), rh = rep(1, Ls),
       model = c("none", "loglinear", "linear", "random"),
       model.control = list())
}
\arguments{
  \item{x}{ A vector containing the values of the stratification variable \eqn{X} for every unit in the population. }
  \item{n}{ A numeric: the target sample size. It has no default value. The argument \code{n} or the argument \code{CV} must be input. }
  \item{CV}{ A numeric: the target coefficient of variation. It has no default value. The argument \code{CV} or the argument \code{n} must be input. }
  \item{Ls}{ A numeric: the number of sampled strata (take-none and certain strata are not counted in \code{Ls}, but here no take-none stratum can be added to the stratified design so \code{Ls} is in fact always equal to \eqn{L}). The default is 3. }
  \item{certain}{ A vector giving the position, in the vector \code{x}, of the units that must be included in the sample (see \code{\link{stratification-package}}). By default \code{certain} is \code{NULL}, which means that no units are chosen a priori to be in the sample. }
  \item{alloc}{ A list specifying the allocation scheme. The list must contain 3 numerics for the 3 exponents \code{q1}, \code{q2} and \code{q3} in the general allocation scheme (see \code{\link{stratification-package}}). The default is Neyman allocation (\code{q1}=\code{q3}=0.5 and \code{q2}=0)}
  \item{rh}{ A vector giving the anticipated response rates in each of the \code{Ls} sampled strata. A single number can be given if the rates do not vary among strata. The default is 1 in each stratum. }
  \item{model}{ A character string identifying the model used to describe the discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can be \code{"none"} if one assumes \eqn{Y=X}, \code{"loglinear"} for the loglinear model with mortality, \code{"linear"} for the heteroscedastic linear model or \code{"random"} for the random replacement model (see \code{\link{stratification-package}} for a description of these models). The default is \code{"none"}.}
  \item{model.control}{ A list of model parameters (see \code{\link{stratification-package}}). The default values of the parameters correspond to the model \eqn{Y=X}. }
  \item{nclass}{ A numeric for the cumulative root frequency method only: the number of classes (Dalenius and Hodges 1959). The default (see \bold{Details}) is \code{min(Ls*15, Nu)} where \code{Nu} is the number of unique values in the \code{x}-vector from which units in the certainty stratum, if any, heve been removed. }
}
\value{
  \item{bh }{A vector of the \eqn{L-1} stratum boundaries proposed by the method.}
  \item{nclassh}{ A vector for the cumulative root frequency method only: the number of classes in each stratum (Dalenius and Hodges 1959). }
  \item{Nh }{A vector of length \eqn{L} containing the population sizes \eqn{N_h}{Nh}, i.e. the number of units in each stratum.}
  \item{nh }{A vector of length \eqn{L} containing the sample sizes \eqn{n_h}{nh}, i.e. the number of units to sample in each stratum. See \code{\link{stratification-package}} for information about the rounding used to get these integer values.}
  \item{n }{The total sample size (\code{sum(nh)}).}
  \item{nhnonint }{A vector of length \eqn{L} containing the non-integer values of the sample sizes, obtained directly from applying the allocation rule (see \code{\link{stratification-package}}).}
  \item{certain.info }{A vector giving statistics for the certainty stratum (see \code{\link{stratification-package}}). It contains \code{Nc}, the number of units chosen a priori to be in the sample, and \code{meanc}, the anticipated mean of \eqn{Y} for these units. }
  \item{opti.nh }{The final value of the criteria to optimize (either the total sample size \eqn{n} if a target \code{CV} was given or the RRMSE if a target \code{n} was given) calculated with the integer stratum sample sizes \code{nh}.}
  \item{opti.nhnonint }{The final value of the criteria to optimize (either the total sample size \eqn{n} if a target \code{CV} was given or the RRMSE if a target \code{n} was given) calculated with the non-integer stratum sample sizes \code{nhnonint}.}
  \item{meanh }{A vector of length \eqn{L} containing the anticipated means of \eqn{Y} in each stratum.}
  \item{varh }{A vector of length \eqn{L} containing the anticipated variances of \eqn{Y} in each stratum.}
  \item{mean }{A numeric: the anticipated global mean value of \eqn{Y}.}
  \item{stderr }{A numeric: the standard error of the anticipated global mean of \eqn{Y}.}
  \item{CV}{The anticipated coefficient of variation for the mean of \eqn{Y}, i.e. \code{stderr} divided \code{mean}. }
  \item{stratumID}{A factor, having the same length as the input \code{x}, which values are either 1, 2, ..., \eqn{L} or \code{"certain"}. The value \code{"certain"} is given to units a priori chosen to be in the sample. This factor identifies, for each observation, the stratum to which it has been assigned.}  
  \item{takeall }{The number of take-all strata in the final solution. Note: It is possible that \eqn{n_h=N_h} for non take-all strata because the condition for an automatic addition of a take-all stratum is \eqn{n_h>N_h}. }
  \item{call }{The function call (object of class "call").}
  \item{date }{A character string that contains the system date and time when the function ended. }
  \item{args }{A list of all the argument values input to the function or set by default. }
}
\details{
The efficiency of the cumulative root frequency method depends on the number of classes \code{nclass} (see Dalenius and Hodges (1959) for a description of these classes). However, there is no theory about how to choose the best value for \code{nclass} (Hedlin 2000). This is a limit of the method.
}
\references{
Baillargeon, S. and  Rivest L.-P. (2011). The construction of stratified designs in R with the package stratification. \emph{Survey Methodology}, \bold{37}(1), 53-65.

Dalenius, T. and Hodges, J.L., Jr. (1959). Minimum variance stratification. \emph{Journal of the American Statistical Association}, \bold{54}, 88-101.

Gunning, P. and Horgan, J.M. (2004). A new algorithm for the construction of stratum boundaries in skewed populations. \emph{Survey Methodology}, \bold{30}(2), 159-166.

Hedlin, D. (2000). A procedure for stratification by an extended Ekman rule. \emph{Journal of Official Statistics}, \bold{61}, 15-29.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{print.strata}}, \code{\link{plot.strata}}, \code{\link{strata.LH}} }
\examples{
### Example for strata.cumrootf
res <- matrix(NA, nrow=20, ncol=2)
i <- 1
for ( n in seq(100,2000,100)){
    cum <- strata.cumrootf(x=MRTS, CV=0.01, Ls=4, alloc=c(0.5,0,0.5), nclass=n)
    res[i,] <- c(n,cum$n)
    i <- i + 1
}
plot(res, ylab="suggested sample size n", xlab="number of classes", main=expression(
     paste("Example of the effect of nclass on n for the cum",sqrt(f)," method")))

### Example for strata.geo
strata.geo(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.35,0.35,0), model="none")
strata.geo(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.35,0.35,0), model="loglinear",
       model.control=list(beta=1.058355, sig2=0.06593083, ph=1))
strata.geo(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.35,0.35,0), rh=0.85,
       model="loglinear", model.control=list(beta=1.058355, sig2=0.06593083, ph=1))
# When non-response or a model is added, the stratum boundaries do not change, 
# only the nh's do.

### Exemple of how a certainty stratum can be usefull with these methods
strata.cumrootf(x=Sweden$REV84, CV=0.05, Ls=4, alloc=c(0.35,0.35,0), model="none",
                nclass=50)
strata.cumrootf(x=sort(Sweden$REV84), CV=0.05, Ls=4, alloc=c(0.35,0.35,0), 
                certain=282:284, model="none", nclass=50)
# The certainty stratum is used here to ensure that the three large units in the
# Sweden$REV84 population are in the sample, since no take-all stratum can be forced 
# in the stratified design with the cumulative root frequency or geometric method.
# We see that this allows to reduce by more than half the suggested sample size n
# (47 vs 19). This example was presented in Baillargeon and Rivest (2011). 

}
\keyword{ survey }
