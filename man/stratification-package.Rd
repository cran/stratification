\name{stratification-package}
\alias{stratification-package}
\alias{stratification}
\docType{package}
\title{ Collection of Functions for Univariate Stratification of Survey Populations }
\description{
This package contains various functions for univariate stratification of survey populations. The well known cumulative root frequency rule by Dalenius and Hodges (1959) and the geometric rule by Gunning and Horgan (2004) are implemented. However, the main function implements a generalized Lavallee-Hidiroglou (1988) method of strata construction. It can be used with Sethi's (1963) or Kozak's (2004) algorithm. The generalized method takes into account a discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. The method can consider a loglinear model with mortality between the variables (Baillargeon, Rivest and Ferland, 2007). When Kozak's algorithm is used, two additional models are available: a heteroscedastic linear model and a random replacement model as in Rivest (2002). The optimal boundaries determination also incorporates, if desired, an anticipated non-response, a take-all stratum for the large units and a take-none stratum for the small units. Moreover, units can be forced to be part of the sample by specifying a certainty stratum. 
}
\details{
\tabular{ll}{
Package: \tab stratification\cr
Type: \tab Package\cr
Version: \tab 2.2-2\cr
Date: \tab 2012-11-03\cr
License: \tab GPL-2\cr
}

\bold{OVERWIEW OF THE FUNCTIONS}\cr

To determine the stratum sample sizes given a set of stratum boundaries: \code{\link{strata.bh}} 

To determine the stratum boundaries first and, in a second step, the stratum sample sizes:\cr 
\code{\link{strata.cumrootf}}: cumulative root frequency method by Dalenius and Hodges (1959)\cr
\code{\link{strata.geo}}: geometric method by Gunning and Horgan (2004)

To determine the optimal stratum boundaries and sample sizes in a single step:\cr
\code{\link{strata.LH}}: generalized Lavallee-Hidiroglou method with Sethi's (1963) or Kozak's (2004) algorithm\cr

All these functions create an object of class "strata", which can be visualized with the S3 methods \code{\link{print.strata}} and \code{\link{plot.strata}}. One can also apply, with the function \code{\link{var.strata}}, a stratified design to a survey variable \eqn{Y} different from the one used for the construction of the stratified design.\cr


\bold{INFORMATION RELATIVE TO MANY FUNCTIONS}\cr

The functions \code{\link{strata.bh}}, \code{\link{strata.cumrootf}}, \code{\link{strata.geo}} and \code{\link{strata.LH}} need to be given:\cr
\code{x}, the values of the stratification variable \eqn{X},\cr
\code{Ls}, the desired number of sampled strata,\cr
\code{alloc}, an allocation rule, and\cr
a target sample size \code{n} or a target level of precision \code{CV} for the survey estimator.\cr 
However, for Sethi's (1963) algorithm, only a target \code{CV} can be given. To reach a target \code{n} using the generalized Lavallee-Hidiroglou method, Kozak's (2004) algorithm has to be used with the \code{\link{strata.LH}} function.\cr   

TYPE OF STRATUM\cr
In this package, four types of stratum exist: take-some, take-none, take-all and certainty. A take-some stratum is a stratum in which some units are sampled. A take-none stratum is a stratum for the smallest units in which no units are sampled. Its purpose is to ignore very small units. On the other hand, a take-all stratum is a stratum for the largest units in which every units are sampled. It allows to insure that the biggest units are in the sample. The following paragraph explains what the special stratum type called \dQuote{certainty} is.\cr  

DEFINITION OF THE CERTAINTY STRATUM\cr
It is possible to insure that some specific units are included in the sample with the argument \code{certain}. This argument is a vector containing the positions in the vector \code{x} of the units to be included with certainty in the sample. We say that these units form the certainty stratum. They are excluded from the population prior to the determination of the stratum boundaries, but they are accounted for in the calculation of the anticipated mean, the RRMSE, the total sample size and the optimization criteria. Essentially, these units form their own separate take-all stratum that is not subject to stratification. They do not have to be consecutive units according to the stratification variable, therefore their variance is meaningless. Non-response is not possible in the certainty stratum. The functions return a value named \code{certain.info} containing the number of units in the certainty stratum and their anticipated mean.\cr   

NUMBER OF STRATA\cr
The \code{Ls} argument represents to the number of sampled strata. The term \dQuote{sampled strata} refers to take-some and take-all strata only. Therefore, take-none and certain strata are not counted in \code{Ls}. If the stratified design does not have a take-none stratum then \code{Ls}=\eqn{L} is the total number of strata, otherwise \code{Ls}=\eqn{L-1}. In the total number of strata \eqn{L}, the certainty stratum, if any, is not counted since we do not need to find its boundaries.\cr

STRATUM NUMBERING\cr
Throughout the package, strata number 1 contains the smallest units and strata number \eqn{L} the biggest ones. So every vector of boundaries contains numbers in ascending order. The function \code{\link{strata.bh}} must be given boundaries \code{bh} fulfilling this condition. This remark also applies to the argument \code{initbh} of \code{\link{strata.LH}} used to give initial boundaries for the optimization algorithm. If a take-none stratum is requested, it is always the first one. On the other hand, if a take-none stratum is requested, it is always the last one. \cr

DEFINITION OF STRATUM BOUNDARIES\cr
Let's note \eqn{b_0, b_1,\ldots, b_L} the stratum boundaries. Stratum \eqn{h} contains all the
units with an \eqn{X}-value in the interval \eqn{[b_{h-1},b_h)} for \eqn{h=1,\ldots,L} such that \eqn{b_0=min(X)} and \eqn{b_L=max(X)+1}, where
\eqn{min(X)} and \eqn{max(X)} are respectively the minimum and the maximum values of the stratification variable. The argument \code{bh} of \code{strata.bh}, the argument \code{initbh} of \code{strata.LH} and the output value \code{bh} of any function of the package \pkg{stratification} with the prefix "strata" are length \eqn{L-1} vectors of the boundaries \eqn{b_1, b_2,\ldots, b_{L-1}}.\cr

DETAILS ABOUT THE TAKE-NONE STRATUM\cr
A non empty take-none stratum induces a bias in the estimator of the mean of \eqn{Y}, and the precision is measured by the relative root mean squared error (RRMSE), not by the coefficient of variation (CV). Regardless, in the functions the argument given to specify a target precision for the survey estimator is always named \code{CV}. However, in the output, the anticipated level of precision is named RRMSE for the functions accepting a \code{takenone} argument (\code{\link{strata.bh}} and \code{\link{strata.LH}}), and it is named CV for the other functions (\code{\link{strata.cumrootf}} and \code{\link{strata.geo}}).

When a \code{takenone} stratum is requested, one can specify a \code{bias.penalty} argument. We define the mean squared error for the estimator of the mean of \eqn{Y} by \eqn{MSE = (bias.penalty \times bias)^2 + variance}{MSE = (bias.penalty x bias)^2 + variance}. It is sometimes possible to estimate the bias using the sum of the \eqn{Y} values in the take-none stratum from administrative data. In this situation, it might be appropriate to set \code{bias.penalty} to a value lower than 1. This will typically enlarge the take-none stratum. The value given to \code{bias.penalty} depends on the confidence level we have in the bias estimate. By default, it is assumed that no bias estimate is available and the whole bias contributes to the MSE (\code{bias.penalty}=1).\cr

SPECIFICATION OF THE ALLOCATION RULE\cr
The \code{alloc} argument must be a list containing the numeric objects \code{q1}, \code{q2} and \code{q3} which specify the allocation rule according to the general allocation scheme presented in Hidiroglou and Srinath (1993) 
\deqn{a_h = \frac{\gamma_h}{\sum_{\mbox{take-some}} \gamma_h} \qquad \mbox{where} \qquad     
\gamma_h=N_h^{2q_1}\bar{Y}_h^{2q_2}S_{yh}^{2q_3}.}{ah = gammah/sum(gammah)     where     gammah = Nh^(2q1) meanYh^(2q2) varYh^(q3).}
Stratum sample sizes are calculated as :
\deqn{
{n_h}_{\mbox{nonint}} = \left\{ \begin{array}{ll}
0 & \mbox{for take-none strata}\\
n \times a_h & \mbox{for take-some strata}\\
N_h & \mbox{for take-all strata}\end{array} \right.
}{nhnonint = 0 for take-none strata, n*ah for take-some strata, Nh for take-all strata}
A proportional allocation is obtained when \code{q1}=0.5 and \code{q2}=\code{q3}=0,\cr
a power allocation is obtained when \code{q1}=\code{q2}=\eqn{p/2} and \code{q3}=0, and\cr
a Neyman allocation (the default) is obtained when \code{q1}=\code{q3}=0.5 and \code{q2}=0.\cr

ROUNDING of the stratum sample sizes \cr
Applying the allocation rule above gives real (non-integer) values for the sample sizes. These are named \code{nhnonint} in the package. The \code{nhnonint} values have to be rounded to get the integer sample sizes, named \code{nh} in the package. Here is how the rounding is done. If a target \code{CV} is requested, the values are simply rounded to the largest integer. However, if a target \code{n} is requested, the rounding is a little more complicated because the \code{nh} should sum to the target \code{n} and we do not want positive nh inferior to 1 to be rounded to zero. Therefore, we first round to 1 the positive nh inferior to 1. Then we calculate how many values (say \code{nup}) must be rounded to the largest integer and how many must be rounded to the smallest integer in order to fulfill the condition \code{sum(nh)=n}. We choose the \code{nup} values with the largest decimal part for the ceiling rounding, the other \code{nh} are rounded down.\cr    

ADJUSTMENT FOR A TAKE-ALL STRATUM\cr
If, after applying the allocation rule, the stratified design contains at least one take-some stratum with \eqn{{n_h}_{\mbox{nonint}} > N_h}{nhnonint>Nh}, the allocation is done again setting the take-some stratum with the largest units as a take-all stratum. This is done until \eqn{{n_h}_{\mbox{nonint}} \leq N_h}{nhnonint<=Nh} for all the take-some strata or until there is only one take-some stratum left. This adjustment is done automatically throughout the package because the target n or CV might not be reached if one omits to do it. Only the function \code{strata.bh} allows not to do it (argument \code{takeall.adjust}).\cr

Note: In special circumstances, the algorithm might result in more than one take-all stratum. If the non-response rate does not vary among the take-all strata, we can see them as forming one big take-all stratum. Otherwise, their boundaries influence the value of the optimization criteria (\eqn{n} or \eqn{CV}). So in the case of a varying non-response rate among the take-all strata, we cannot see them as forming one big take-all stratum.\cr

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
  \item{\code{pcertain}}{ A numeric: the survival rate in the certainty stratum, if a certainty stratum is added to the stratified design. The default is 1. }
  \item{\code{gamma}}{ A numeric: the exponent of \eqn{X} in the residual variance of the "linear" model. The default is 0. }
  \item{\code{epsilon}}{ A numeric: the probability that the \eqn{Y}-value for a unit is equal to the \eqn{X}-value for a randomly selected unit in the population. It concerns the "random" model only. The default is 0.\cr}
}
Note: The default values of the parameters simplify any model to \eqn{Y=X}. Therefore, the default is always to consider that there is no discrepancy between the stratification and the survey variables. The \code{model} argument even has the default value \code{"none"}, which also means \eqn{Y=X}.

}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\references{
Baillargeon, S., Rivest, L.-P., Ferland, M. (2007). Stratification en enquetes entreprises : Une revue et quelques avancees. \emph{Proceedings of the Survey Methods Section, 2007 SSC Annual Meeting}. \url{ http://www.ssc.ca/survey/documents/SSC2007_S_Baillargeon.pdf }

Baillargeon, S. and Rivest, L.-P. (2009). A general algorithm for univariate stratification. \emph{International Stratification Review}, \bold{77}(3), 331-344.

Baillargeon, S. and  Rivest L.-P. (2011). The construction of stratified designs in R with the package stratification. \emph{Survey Methodology}, \bold{37}(1), 53-65. \url{ http://www.statcan.gc.ca/pub/12-001-x/2011001/article/11447-eng.pdf }

Dalenius, T. and Hodges, J.L., Jr. (1959). Minimum variance stratification. \emph{Journal of the American Statistical Association}, \bold{54}, 88-101.

Gunning, P. and Horgan, J.M. (2004). A new algorithm for the construction of stratum boundaries in skewed populations. \emph{Survey Methodology}, \bold{30}(2), 159-166.

Hidiroglou, M.A. and Srinath, K.P. (1993). Problems associated with designing subannual business surveys. \emph{Journal of Business & Economic Statistics}, \bold{11}, 397-405.

Kozak, M. (2004). Optimal stratification using random search method in agricultural surveys. \emph{Statistics in Transition}, \bold{6}(5), 797-806.

Lavallee, P. and Hidiroglou, M.A. (1988). On the stratification of skewed populations. \emph{Survey Methodology}, \bold{14}, 33-43.

Rivest, L.-P. (2002). A generalization of the Lavallee and Hidiroglou algorithm for stratification in business surveys. \emph{Survey Methodology}, \bold{28}(2), 191-198.

Sethi, V. K. (1963). A note on optimum stratification of populations for estimating the population means. \emph{The Australian Journal of Statistics}, \bold{5}, 20-33.
}
\keyword{ package }
