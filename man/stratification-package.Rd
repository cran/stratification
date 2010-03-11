\name{stratification-package}
\alias{stratification-package}
\alias{stratification}
\docType{package}
\title{ Collection of Functions for Univariate Stratification of Survey Populations }
\description{
This package contains various functions for univariate stratification of survey populations. The well known cumulative root frequency rule by Dalenius and Hodges (1959) and the geometric rule by Gunning and Horgan (2004) are implemented. However, the main function implements a generalized Lavallee-Hidiroglou (1988) method of strata construction. It can be used with Sethi's (1963) or Kozak's (2004) algorithm. The generalized method takes into account a discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. The method can consider a loglinear model with mortality between the variables. When Kozak's algorithm is used, two additional models are available: a heteroscedastic linear model and a random replacement model as in Rivest (2002). The optimal boundaries determination also incorporates, if desired, an anticipated non-response, a take-all stratum for the large units and a take-none stratum for the small units.
}
\details{
\tabular{ll}{
Package: \tab stratification\cr
Type: \tab Package\cr
Version: \tab 2.0-0\cr
Date: \tab 2010-03-11\cr
License: \tab GPL-2\cr
}

\bold{OVERWIEW OF THE FUNCTIONS}\cr

To determine the stratum sample sizes given a set of stratum boundaries: \code{\link{strata.bh}} 

To determine the stratum boundaries first and, in a second step, the stratum sample sizes:\cr 
\code{\link{strata.cumrootf}}: cumulative root frequency method by Dalenius and Hodges (1959)\cr
\code{\link{strata.geo}}: geometric method by Gunning and Horgan (2004)

To determine the optimal stratum boundaries and sample sizes in a single step:\cr
\code{\link{strata.LH}}: generalized Lavallee-Hidiroglou method with Sethi's (1963) or Kozak's (2004) algorithm\cr

All these functions create objects of class "strata" which can be visualized with \code{\link{print.strata}} and \code{\link{plot.strata}}. One can also apply, with the function \code{\link{var.strata}}, a stratified design to a survey variable \eqn{Y} different from the one used for the construction of the stratified design.\cr


\bold{INFORMATION RELATIVE TO MANY FUNCTIONS}\cr

The functions \code{\link{strata.bh}}, \code{\link{strata.cumrootf}}, \code{\link{strata.geo}} and \code{\link{strata.LH}} need to be given \code{x}, the values of the stratification variable \eqn{X}, \code{Ls}, the desired number of strata, \code{alloc}, an allocation rule, and a target sample size \code{n} or a target level of precision \code{CV} for the survey estimator. However, for Sethi's (1963) algorithm, only a target \code{CV} can be given. To reach a target \code{n} using the generalized Lavallee-Hidiroglou method, Kozak's (2004) algorithm has to be used with the \code{\link{strata.LH}} function.\cr   

Throughout the package, strata number 1 contains the smallest units and strata number \eqn{L} the biggest ones. So every vector of boundaries contains numbers in ascending order. The function \code{\link{strata.bh}} must be given boundaries \code{bh} fulfilling this condition. This remark also applies to the argument \code{initbh} of \code{\link{strata.LH}} used to give initial boundaries.\cr

MORE DETAILS ABOUT THE TAKE-NONE STRATUM\cr
The \code{Ls} argument refers to the number of sampled strata. If the stratified design does not have a take-none stratum then \code{Ls}=\eqn{L} is the total number of strata, otherwise \code{Ls}=\eqn{L-1}.

A non empty take-none stratum induces a bias in the estimator of the mean of \eqn{Y}, and the precision is measured by the relative root mean squared error (RRMSE), not by the coefficient of variation (CV). Regardless, in the functions the argument given to specify a target precision for the survey estimator is always named \code{CV}. However, in the output, the anticipated level of precision is named RRMSE for the functions accepting a \code{takenone} argument (\code{\link{strata.bh}} and \code{\link{strata.LH}}), and it is named CV for the other functions (\code{\link{strata.cumrootf}} and \code{\link{strata.geo}}).

When a \code{takenone} stratum is requested, one can specify a \code{bias.penalty} argument. We define the mean squared error for the estimator of the mean of \eqn{Y} by \eqn{MSE = (bias.penalty \times bias)^2 + variance}{MSE = (bias.penalty x bias)^2 + variance}. It is sometimes possible to estimate the bias using the sum of the \eqn{Y} values in the take-none stratum from administrative data. In this situation, it might be appropriate to set \code{bias.penalty} to a value lower than 1. This will typically enlarge the take-none stratum. The value given to \code{bias.penalty} depends on the confidence level we have in the bias estimate. By default, it is assumed that no bias estimate is available and the whole bias contributes to the MSE (\code{bias.penalty}=1).\cr


CERTAIN STRATUM\cr
It is possible to insure that some units are included in the sample with the argument \code{certain}. This argument is a vector containing the positions in the vector \code{x} of the units to be included with certainty in the sample. We say that these units form the certain stratum. They are excluded from the population prior to the determination of the stratum boundaries, but they are accounted for in the calculation of the anticipated mean and the RRMSE. Essentially, these units form their own separate take-all stratum that is not subject to stratification. They do not have to be consecutive units according to the stratification variable, therefore their variance is meaningless. Non-response is not possible in the certain stratum. The functions return a value named \code{certain.info} containing the number of units in the certain stratum and their anticipated mean.\cr   

SPECIFICATION OF THE ALLOCATION RULE\cr
The \code{alloc} argument must be a list containing the numerics \code{q1}, \code{q2} and \code{q3} which specify the allocation rule according to the general allocation scheme presented in Hidiroglou and Srinath (1993) \deqn{a_h = \frac{\gamma_h}{\sum \gamma_h} \qquad \mbox{where} \qquad     \gamma_h=N_h^{2q_1}\bar{Y}_h^{2q_2}S_{yh}^{2q_3}.}{ah = gammah/sum(gammah)     where     gammah = Nh^(2q1) meanYh^(2q2) varYh^(q3).} For example, a proportional allocation is obtained when \code{q1}=0.5 and \code{q2}=\code{q3}=0, a power allocation is obtained when \code{q1}=\code{q2}=\eqn{p/2} and \code{q3}=0, and a Neyman allocation (the default) is obtained when \code{q1}=\code{q3}=0.5 and \code{q2}=0.\cr

SPECIFICATION OF A MODEL BETWEEN \eqn{Y} AND \eqn{X}\cr
Every function can take into account a discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. The functions \code{\link{strata.bh}}, \code{\link{strata.cumrootf}} and \code{\link{strata.geo}} perform allocation on the basis of anticipated moments whereas the \code{\link{strata.LH}} function goes further; it determines the optimal boundaries considering the anticipated moments. The following models for the relationship between \eqn{Y} and \eqn{X} can be specified through the \code{model} and \code{model.control} arguments:\cr

\bold{- loglinear model with mortality} (\code{model="loglinear"}):
\deqn{
 Y = \left\{ \begin{array}{ll}
\exp(\alpha + \mbox{\code{beta}} \ \log(X) + \mbox{\code{epsilon}}) & \mbox{with probability } p_h\\
0 & \mbox{with probability } 1-p_h \end{array} \right.
}{Y = exp(alpha + beta log(X) + epsilon) with probability ph, 0 with probability 1-ph}
where \eqn{\mbox{\code{epsilon}} \sim N(0,\mbox{\code{sig2}})}{epsilon ~ N(0,sig2)} is independent of \eqn{X}. The parameter \eqn{p_h}{ph}
is specified through \code{ph}, \code{ptakenone} and \code{pcertain} (elements of \code{model.control}). Note: The \eqn{\alpha}{alpha} parameter does not have to be specified because \eqn{exp(\alpha)}{exp(alpha)} is a multiplicative factor that has no impact on the outcome.\cr

\bold{- heteroscedastic linear model} (\code{model="linear"}):
\deqn{
Y = \mbox{\code{beta}} X + \mbox{\code{epsilon}}
}{ Y = beta X + epsilon}
where \eqn{\mbox{\code{epsilon}} \sim N(0,\mbox{\code{sig2}} \ X^{\mbox{\code{gamma}}})}{epsilon ~ N(0,sig2 X^gamma)}.\cr

\bold{- random replacement model} (\code{model="random"}):
\deqn{
 Y = \left\{ \begin{array}{ll}
X & \mbox{with probability } 1-\mbox{\code{epsilon}} \\
Xnew & \mbox{with probability } \mbox{\code{epsilon}} \end{array} \right.
}{Y = X with probability 1-epsilon, Xnew with probability epsilon}
where \eqn{Xnew} is a random variable independent of \eqn{X} having the same distribution than \eqn{X}.\cr

The \code{model.control} argument is a list that can supply any of the following model parameter:
\describe{
  \item{\code{beta}}{ A numeric: the slope of the "loglinear" or "linear" model. The default is 1. }
  \item{\code{sig2}}{ A numeric: the variance parameter of the "loglinear" or "linear" model. The default is 0. }
  \item{\code{ph}}{ A vector giving the survival rate in each of the \code{Ls} sampled strata for the "loglinear" model. A single number can be given if the rate doesn't vary among strata. The default is 1 in each stratum. }
  \item{\code{ptakenone}}{ A numeric: the survival rate in the take-none stratum, if a take-none stratum is added to the stratified design. The default is 1. }
  \item{\code{pcertain}}{ A numeric: the survival rate in the certain stratum, if a certain stratum is added to the stratified design. The default is 1. }
  \item{\code{gamma}}{ A numeric: the exponent of \eqn{X} in the residual variance of the "linear" model. The default is 0. }
  \item{\code{epsilon}}{ A numeric: the probability that the \eqn{Y}-value for a unit is equal to the \eqn{X}-value for a randomly selected unit in the population. It concerns the "random" model only. The default is 0.}
}

.\cr

ADJUSTMENT FOR A TAKE-ALL STRATUM\cr
If, at the end of the computations, the solution reached by any of these functions contains at least one take-some stratum with \eqn{n_h > N_h}{nh>Nh}, the computations are carried out again setting the take-some stratum with the largest units as a take-all stratum. This is done until \eqn{n_h \leq N_h}{nh<=Nh} for all the take-some strata. For the function \code{strata.bh} only, this adjustment is optional (argument \code{takeall.adjust}).\cr

ROUNDING of the \eqn{n_h}{nh}'s \cr
These functions return the sample sizes \code{nh}, which are integers, but also the non-integer values of these sample sizes, named \code{nh.nonint}. The vector \code{nh.nonint} is the result of the application of the allocation rule. The \code{nh.nonint} values are rounded to get the integer sample sizes \code{nh}. Here is how the rounding is done. If a target \code{CV} is requested, the values are simply rounded to the largest integer. However, if a target \code{n} is requested, the rounding is a little more complicated because the \code{nh} should sum to the target \code{n} and we do not want positive nh inferior to 1 to be rounded to zero. Therefore, we first round to 1 the positive nh inferior to 1. Then we calculate how many values (say \code{nup}) must be rounded to the largest integer and how many must be rounded to the smallest integer in order to fulfill the condition \code{sum(nh)=n}. We choose the \code{nup} values with the largest decimal part for the ceiling rounding, the other \code{nh} are rounded down.    

}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\references{
Baillargeon, S. and Rivest, L.-P. (2010). Univariate stratification in R with the package stratification. Paper in progress.

Baillargeon, S., Rivest, L.-P., Ferland, M. (2007). Stratification en enquetes entreprises : Une revue et quelques avancees. \emph{Proceedings of the Survey Methods Section, 2007 SSC Annual Meeting}. \url{ http://www.ssc.ca/survey/documents/SSC2007_S_Baillargeon.pdf }

Dalenius, T. and Hodges, J.L.,Jr. (1959). Minimum variance stratification. \emph{Journal of the American Statistical Association}, \bold{54}, 88-101.

Gunning, P. and Horgan, J.M. (2004). A new algorithm for the construction of stratum boundaries in skewed populations. \emph{Survey Methodology}, \bold{30}(2), 159-166.

Hidiroglou, M.A. and Srinath, K.P. (1993). Problems associated with designing subannual business surveys. \emph{Journal of Business & Economic Statistics}, \bold{11}, 397-405.

Kozak, M. (2004). Optimal stratification using random search method in agricultural surveys. \emph{Statistics in Transition}, \bold{6}(5), 797-806.

Lavallee, P. and Hidiroglou, M.A. (1988). On the stratification of skewed populations. \emph{Survey Methodology}, \bold{14}, 33-43.

Rivest, L.-P. (2002). A generalization of the Lavallee and Hidiroglou algorithm for stratification in business surveys. \emph{Survey Methodology}, \bold{28}(2), 191-198.

Rivest, L.-P. and Baillargeon, S. (2009). A general algorithm for univariate stratification. \emph{International Stratification Review}, \bold{77}(3), 331-344.

Sethi, V. K. (1963). A note on optimum stratification of populations for estimating the population means. \emph{The Australian Journal of Statistics}, \bold{5}, 20-33.
}
\keyword{ package }
