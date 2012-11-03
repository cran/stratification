\name{strata.LH}
\alias{strata.LH}
\title{ Generalized Lavallee-Hidiroglou Method of Strata Construction }
\description{
This function aims at constructing optimal strata with a generalized Lavallee-Hidiroglou (1998) method.
The function uses Kozak's (2004) algorithm by default, but it can also apply Sethi's (1963) algorithm (argument \code{algo="Sethi"}). The function can take into account a discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can consider a loglinear model with mortality between the variables (Baillargeon and Rivest, 2009). With Kozak's algorithm, two additional models are implemented: an heteroscedastic linear model and a random replacement model as in Rivest (2002). The determination of the optimal boundaries also incorporates, if desired, an anticipated non-response, a take-all stratum for the large units, a take-none stratum for the small units, and a certainty stratum to ensure that some specific units are in the sample.\cr\cr
Sethi's algorithm is not used by default because it can be numerically unstable, especially with a take-none stratum. Better results were obtained with Kozak's algorithm in our numerical experiments.
}
\usage{
strata.LH(x, n = NULL, CV = NULL, Ls = 3, certain = NULL,
          alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), takenone = 0, 
          bias.penalty = 1, takeall = 0, rh = rep(1, Ls), 
          model = c("none", "loglinear", "linear", "random"), 
          model.control = list(), initbh = NULL, 
          algo = c("Kozak", "Sethi"), algo.control = list())
}
\arguments{
  \item{x}{ A vector containing the values of the stratification variable \eqn{X} for every unit in the population. }
  \item{n}{ A numeric: the target sample size. It has no default value. The argument \code{n} or the argument \code{CV} must be input. }
  \item{CV}{ A numeric: the target coefficient of variation or relative root mean squared error if \code{takenone}=1. It has no default value. The argument \code{CV} or the argument \code{n} must be input. }
  \item{Ls}{ A numeric: the number of sampled strata (take-none and certain strata are not counted in \code{Ls}). The default is 3. }
  \item{certain}{ A vector giving the position, in the vector \code{x}, of the units that must be included in the sample (see \code{\link{stratification-package}}). By default \code{certain} is \code{NULL}, which means that no units are chosen a priori to be in the sample. }
  \item{alloc}{ A list specifying the allocation scheme. The list must contain 3 numerics for the 3 exponents \code{q1}, \code{q2} and \code{q3} in the general allocation scheme (see \code{\link{stratification-package}}). The default is Neyman allocation (\code{q1}=\code{q3}=0.5 and \code{q2}=0)}
  \item{takenone}{ A numeric: the number of take-none strata (0 or 1). The default is 0, i.e. no take-none stratum is included.}
  \item{bias.penalty}{ A numeric between 0 and 1 giving the penalty for the bias in the anticipated mean squared error (MSE) of the survey estimator (see \code{\link{stratification-package}}). This argument is relevant only if \code{takenone}=1. The default is 1.}
  \item{takeall}{ A numeric: the number of take-all strata (one of \{0, 1, ..., \code{Ls}-1\}). The default is 0, i.e. no take-all stratum is included.}
  \item{rh}{ A vector giving the anticipated response rates in each of the \code{Ls} sampled strata. A single number can be given if the rates do not vary among strata. The default is 1 in each stratum. }
  \item{model}{ A character string identifying the model used to describe the discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can be \code{"none"} if one assumes \eqn{Y=X}, \code{"loglinear"} for the loglinear model with mortality, \code{"linear"} for the heteroscedastic linear model or \code{"random"} for the random replacement model (see \code{\link{stratification-package}} for a description of these models). The default is \code{"none"}, so the original Lavallee-Hidiroglou (1998) method of strata construction is used. The last two models \code{"linear"} and \code{"random"} are only available with Kozak's algorithm.}
  \item{model.control}{ A list of model parameters (see \code{\link{stratification-package}}). The default values of the parameters correspond to the model \eqn{Y=X}. }
  \item{initbh}{ A vector of initial stratum boundaries  (see \bold{Details}).}  
  \item{algo}{ A character string identifying which optimization algorithm is to be used. It can be "Kozak" for Kozak's (2004) algorithm or "Sethi" for Sethi's (1963) algorithm (see \bold{Details}). The default is "Kozak" because it performs better and offers more options than Sethi's algorithm. }
  \item{algo.control}{ A list of parameters to control the optimization algorithm (see \bold{Details}).}
}
\details{
CALCULATION OF THE OPTIMIZATION CRITERIA

Kozak's and Sethi's algorithms aim at finding optimal stratum boundaries. This optimality refers to a criteria :\cr
the total sample size \eqn{n} if a target \code{CV} was given\cr
the RRMSE if a target \code{n} was given (not available for Sethi's algorithm).\cr
This criteria, no matter if it is \eqn{n} or the RRMSE, is a function of the stratum sample sizes, among others. But these stratum sample sizes can be non-integer (\code{nhnonint} obtained directly from applying the allocation rule) or integer (\code{nh} obtained from rounding the \code{nhnonint}), as mentioned in \code{\link{stratification-package}}. The optimization criteria does not vary much depending on the sample size used, but it is not exactly the same. Let's call \code{opti.nh} the criteria calculated with \code{nh} and \code{opti.nhnonint} the criteria calculated with \code{nhnonint}. Optimizing based on \code{opti.nh} vs \code{opti.nhnonint} might not give the same result. In fact, according to our numerical tests, it is a little better to optimize based on \code{opti.nh}. This is logical since the sample sizes used in practice, to apply a stratified design, are the integer one. Therefore, as of version 2.1-0 of the package, Kozak's algorithm uses by default the integer sample sizes \code{nh} to calculate the optimization criteria (but it can also use \code{opti.nhnonint} by setting  \code{idopti="nhnonint"} in \code{algo.control}). However, Sethi's algorithm works with derivatives to perform optimization. Therefore, it does not treat the optimization as a discrete problem as Kozak's algorithm does. It can only work with a real optimization criteria, i.e. the criteria \code{opti.nhnonint} calculated with non-integer sample size. This is one more reason to favour Kozak's algorithm over Sethi's one.\cr\cr


INITIAL STRATUM BOUNDARIES

Both Kozak's and Sethi's algorithm are iterative and they need starting values. It can be given by the user with the \code{initbh} argument. If not given, the default initial boundaries are obtained with the function \code{\link{strata.cumrootf}} (cumulative root frequency method) for Kozak's algorithm. For Sethi's algorithm, the default value of \code{initbh} is the set of \code{Ls}-1 equidistant points along the range of the \eqn{X}-values (arithmetic starting points of Gunning and Horgan (2007)). Let's define \eqn{d=(max(X)-min(X))/Ls}, the arithmetic boundaries are \eqn{b_h = min(X) + h \times d}{bh = min(X) + h x d} for \eqn{h=1,\ldots,Ls-1}. 

The default initial boundaries are not the same for the two algorithms because Kozak's algorithm uses by default the criteria \code{opti.nh} whereas Sethi's algorithm can only use \code{opti.nhnonint}. In our numerical experiments with the criteria \code{opti.nh}, the cumulative root frequency boundaries performed a little better than geometric, arithmetic or robust boundaries (for a definition of theses boundaries, see details about the \code{trymany} algorithm parameter). On the other hand, the arithmetic boundaries had good performances (Baillargeon and Rivest 2009) with the criteria \code{opti.nhnonint}. Let's note however that we never found, in any of our numerical experiments, initial boundaries clearly better than others.

The length of the vector \code{initbh} can be \code{Ls}-1 or \eqn{L}-1, where \eqn{L} is the total number of strata. When \code{takenone=0}, \eqn{L}=\code{Ls} and it does not make any difference. But for a stratified design with a take-none stratum, it means that one can give as \code{initbh} argument the vector of initial boundaries \eqn{(b_1, b_2, \ldots, b_{L-1})}{(b1, b2, ..., bL-1)} or \eqn{(b_2, \ldots, b_{L-1})}{(b2, ..., bL-1)}. In the second option, the upper boundary of the take-none stratum is not given. In that case, it is set by default to 
the first percentile of \code{x}. If this first percentile is equal to the minimum value of \code{x}, this initial 
boundary would lead to an empty take-none stratum. In that case, the initial boundary of the take-none stratum is 
rather set to the second smallest value of \code{x}.

Sometimes, the specified initial boundaries are not suitable for Kozak's algorithm. As written below, this algorithm 
verifies at each iteration that the sampled strata contains at least \code{minNh} units and that they have positive 
sample sizes \code{nh}. If the initial boundaries do not meet these conditions, the initial optimization criteria is not comparable to the criteria for boundaries respecting the conditions. Therefore, in such a situation, a warning is produced and the algorithm is not run. \cr\cr  


ALGORITHMS

\bold{Sethi:} The formulas implemented for Sethi's algorithm are presented in Baillargeon, Rivest and Ferland (2007). As mentioned previously, this algorithm is available for a target CV only. Moreover, it works with derivatives to perform optimization. Therefore, it can only work with a real optimization criteria (\code{opti.nhnonint}).

\bold{Kozak:} Kozak's algorithm is described in Kozak (2004). It starts at the initial boundaries \code{initbh} and at each iteration, it chooses a stratum boundary at random and a random modification for this boundary among the 2*\code{maxstep} possible alternatives. If this modification reduces the optimization criteria, creates sampled strata containing at least \code{minNh} units and leads to positive \code{nh}, it is accepted. Otherwise, the boundaries are not changed. If the boundaries remain unchanged for \code{maxstill} consecutive iterations, the algorithm has reached convergence and it stops.\cr
Kozak's (2004) used the condition \code{nh} >= 2 instead of \code{nh} > 0. We chose to modify this requirement because we noticed that more severe conditions sometimes prevented the algorithm from selecting a path leading to the optimal solution. According to our numerical experiments, if one sets a posteriori to 2 the \code{nh}'s that are equal to 1 at the end of the algorithm, the new sample size is smaller than or equal to the one obtained with the condition \code{nh} >= 2.\cr
If no \code{initbh} argument is given, the algorithm is run with three different sets of initial boundaries (see details of the algorithm parameter \code{trymany}). If \code{initbh} is given, it uses these initial boundaries. Bot no matter how many initial boundaries are tried, the algorithm is run \code{rep} times with each. Running the algorithm many times helps having more stable results. Since the algorithm is random, two runs of the algorithm, event from the same starting point, can lead to different results. Setting \code{rep} to a large value and trying more than one set of initial boundaries increases the calculation time, but also increases the chance of finding the global minimum.  also helps in finding the global minimum. 

The optimization of stratum boundaries given a target CV or a target \eqn{n} is a discrete problem. The number of possible sets of boundaries is \code{choose(Nu-1,L-1)} when \code{takenone=0}, and \code{choose(Nu,L-1)} when \code{takenone=1}, where \code{Nu} stands for the number of unique values in the \code{x}-vector from which units in the certainty stratum, if any, heve been removed. When this number is not too large, trying them all is numerically possible. That's what the function \code{strata.LH} do, for Kozak's algorithm, when the number of possible solutions is lower than the \code{minsol} parameter in \code{algo.control}. This complete enumeration of all possible cases ensure that the global minimum is reached. 

Note: Since version 2.2-0 of the package \pkg{stratification}, the modified Kozak's algorithm, which is a non-random version of Kozak's algorithm, is no more available because it never performed as well as the original Kozak's algorithm in our numerical experiments.\cr\cr


ALGORITHM PARAMETERS

The \code{algo.control} argument is a list to supply any of the following parameters which control the algorithm. Sethi's algorithm only uses the first argument \code{maxiter}.
\describe{
  \item{\code{maxiter}}{ A numeric: the maximal number of iterations. The default is 500 for Sethi's algorithm and 10 000 for Kozak's. It is only used to prevent from infinite loops in case of non-convergence.  }
  \item{\code{minsol}}{ A numeric for Kozak's algorithm only: if the number of possible sets of boundaries is lower than \code{minsol}, a complete enumeration of the solutions is performed instead of running the algorithm. Every set of boundaries is tried, which ensures that the global minimum is reached. The default value of \code{minsol} is 10 000. This parameter has to take a value between 2 and 2 000 000. }
  \item{\code{idopti}}{ A character string for Kozak's algorithm only: this argument determines which stratum sample sizes are going to be used to calculate the optimization criteria. It can take the value \code{"nh"} (criteria calculated with the integer sample sizes \code{nh}) or \code{"nhnonint"} (criteria calculated with the non-integer sample sizes \code{nhnonint}). The default value is \code{"nh"} since it gives slightly better results than \code{idopti="nhnonint"} and also because the integer sample sizes are the ones used in practice. When a complete enumeration is performed, the criteria is automatically calculated with integer stratum sample sizes. Note: Prior to version 2.1-0 of the package \pkg{stratification}, only the option \code{idopti="nhnonint"} was available (in fact the algorithm parameter \code{idopti} did not exist). }
  \item{\code{minNh}}{ A numeric for Kozak's algorithm only: the minimum number of units required in each sampled stratum (no restriction is put on the take-none stratum, if included). \code{minNh} must be greater or equal to 2, which is the default. }
  \item{\code{maxstep}}{ A numeric for Kozak's algorithm only: the maximal step for boundary modification (see the algorithm description above). The default is \code{Nu/10}, rounded up and truncated to 100 (\code{Nu} is the number of unique values in the \code{x}-vector from which units in the certainty stratum, if any, heve been removed). This default value is, for most populations, much larger than the initial suggestion by Kozak (2004) of a integer no bigger than 5. Our numerical experiments showed us that a \code{maxstep} value of about 3 (package \pkg{stratification} initial default value) is adequate when the optimization criteria is calculated with non-integer stratum sample sizes (\code{opti.nhnonint}). However, when this criteria is calculated with integer stratum sample sizes (\code{opti.nh}), the algorithm never do worst and sometimes reaches an even lower optimization criteria with a larger \code{maxstep}. The parameter \code{maxstill} must be increased consequently (see below).  The downside of large \code{maxstep} and \code{maxstill} values is that many iterations are needed for the algorithm to converge, making it longer to run. We implemented a solution to this problem. This solution is based on the fact that when the algorithm is close to the solution, it can only accept boundaries modifications with small steps. The large steps are useful during the first iterations of the algorithm to let the boundaries move freely, but it becomes useless when convergence is close. Therefore, the values of \code{maxstep} and \code{maxstill} are brought down to 3 and 50, respectively, when new boundaries are accepted under these circumstances : the step in the boundary modification is lower than 3 in absolute value, more than 50 iterations without change happened before this change, the relative change in the optimization criteria is small. This rule is arbitrary, but it seems to work well so far. It gives results as good as without changing the values of \code{maxstep} and \code{maxstill}, but it runs much faster because it needs less iterations to converge.}
  \item{\code{maxstill}}{ A numeric for Kozak's algorithm only: the maximal number of iterations without a change in the boundaries (see the algorithm description above). The default is \code{maxstep}*10, bounded between 50 and 500. It depends on the value of \code{maxstep} because we increased the default value of this parameter by making it a function of the number of unique values in the \code{x}-vector (see above). The upper bound of 500 ensures that the algorithm does not have to make too much iterations (which take time) to converge. }
  \item{\code{rep}}{ A numeric for Kozak's algorithm only: the number of repetitions of the algorithm (see the algorithm description above). The default is 5. Since version 2.2-0 of the package \pkg{stratification}, the option \code{rep="change"} has been replaced by the option \code{trymany=TRUE}.}
  \item{\code{trymany}}{ A logical for Kozak's algorithm only: if \code{trymany=TRUE}, three sets of initial boundaries are tried instead of one. These boundaries are cumulative root frequency (Dalenius and Hodges, 1959), geometric (Gunning and Horgan, 2004) and robust ones. Robust boundaries were created in order to have boundaries respecting as often as possible the conditions imposed on boundaries in Kozak's algorithm: strata containing at least \code{minNh} units and positive \code{nh}. Robust boundaries give an empty take-none stratum if such a stratum is requested, take-all strata as small as possible, and take-some strata with approximately the same number of unique \eqn{X}-values. When \code{trymany=TRUE}, the same values for the other algorithm parameters (the one given by the user or otherwise the default) are used for every set of initial boundaries. If \code{trymany=FALSE}, only the given or default (see above) initial boundaries are used. If a user give an \code{initbh} argument, but also select \code{trymany=TRUE}, the later is ignored and only the \code{initbh} initial boundaries are used.\cr\cr}  
}

CONVERGENCE

It is possible for the algorithm not to converge. In this case, a warning is printed. The only possible cause of non-convergence for Kozak's algorithm is to reach the maximum number of iterations before the stopping rule has been met. For the algorithm to converge, the argument \code{maxiter} has to be increased. On the other hand, non-convergence of Sethi's algorithm has several possible causes:
\itemize{
\item a division by zero caused by an empty stratum can occur ;
\item a division by zero caused by a 0 stratum variance can occur ;
\item the square root of a negative number (negative discriminant) can occur ;
\item the maximum number of iterations can be reached (often because the algorithm is caught in a loop of non-optimal sets of boundaries).
}
If a non-convergence happens, the user can try to change the initial boundaries or the model parameters. The user can also choose to work with Kozak's algorithm which should converge given an appropriate maximum number of iterations.\cr\cr

WARNING ABOUT LOCAL MINIMA

Let's note that even if the algorithm converges, it is not guaranteed that it has reached a global minimum. Local minimum can occur both with Sethi's and Kozak's algorithms (Rivest and Baillargeon, 2009).
}
\value{
  \item{bh }{A vector of the \eqn{L-1} optimal stratum boundaries found by the algorithm.}
  
  \item{Nh }{A vector of length \eqn{L} containing the population sizes \eqn{N_h}{Nh}, i.e. the number of units in each stratum.}
  
  \item{nh }{A vector of length \eqn{L} containing the sample sizes \eqn{n_h}{nh}, i.e. the number of units to sample in each stratum. See \code{\link{stratification-package}} for information about the rounding used to get these integer values.}
  
  \item{n }{The total sample size (\code{sum(nh)}).}
  
  \item{nhnonint }{A vector of length \eqn{L} containing the non-integer values of the sample sizes, obtained directly from applying the allocation rule (see \code{\link{stratification-package}}).}
  
  \item{certain.info }{A vector giving statistics for the certainty stratum (see \code{\link{stratification-package}}). It contains \code{Nc}, the number of units chosen a priori to be in the sample, and \code{meanc}, the anticipated mean of \eqn{Y} for these units. }
  
  \item{opti.criteria }{The final value of the criteria to optimize : either \eqn{n} if a target \code{CV} was given or the RRMSE if a target \code{n} was given. If the algorithm parameter \code{idopti} takes the value \code{"nh"} (the default), the stratum sample sizes used for the calculation of this criteria are the integer ones (\code{nh}). If the algorithm parameter \code{idopti} takes the value \code{"nhnonint"}, the non-integer stratum sample sizes (\code{nhnonint}) are used. In other words, \code{idopti} determines here if \code{opti.criteria} is \code{opti.nh} or \code{opti.nhnonint}.}
  
  \item{meanh }{A vector of length \eqn{L} containing the anticipated means of \eqn{Y} in each stratum.}
  
  \item{varh }{A vector of length \eqn{L} containing the anticipated variances of \eqn{Y} in each stratum.}
  
  \item{mean }{A numeric: the anticipated global mean value of \eqn{Y}.}
  
  \item{RMSE }{A numeric: the root mean squared error (or standard error if \code{takenone}=0) of the anticipated global mean of \eqn{Y}. This is defined as the squared root of: (\code{bias.penalty} x bias of the mean)^2 + variance of the mean.}
  
  \item{RRMSE }{A numeric: the anticipated relative root mean squared error (or coefficient of variation if \code{takenone}=0) for the mean of \eqn{Y}, i.e. \code{RMSE} divided by \code{mean}. }
  
  \item{relativebias }{A numeric: the anticipated relative bias of the estimator, i.e. (\code{bias.penalty} x bias of the mean) divided by \code{mean}. If \code{takenone}=0, this numeric is zero. }
  
  \item{propbiasMSE }{A numeric: the proportion of the MSE attributable to the bias of the estimator, i.e. (\code{bias.penalty} x bias of the mean)^2 divided by the MSE of the \code{mean}. If \code{takenone}=0, this numeric is zero.}
  
  \item{stratumID}{A factor, having the same length as the input \code{x}, which values are either 1, 2, ..., \eqn{L} or \code{"certain"}. The value \code{"certain"} is given to units a priori chosen to be in the sample. This factor identifies, for each observation, the stratum to which it has been assigned.}  
  
  \item{takeall }{The number of take-all strata in the final solution. Note: It is possible that \eqn{n_h=N_h} for non take-all strata because the condition for an automatic addition of a take-all stratum is \eqn{n_h>N_h}. }
  
  \item{call }{The function call (object of class "call").}
  
  \item{date }{A character string that contains the system date and time when the function ended. }
  
  \item{args }{A list of all the argument values input to the function or set by default. }
  
  \item{iter.detail }{Only if an iterative algorithm is run (thus not int the output when a complete enumeration is performed): a data frame giving, for each iteration of the algorithm:\cr
\code{"bh"}: the stratum boundaries;\cr
\code{"opti.nh"}: the value of the criteria to optimize calculated with integer stratum sample sizes (for Kozak's algorithm only);\cr
\code{"opti.nhnonint"}: the value of the criteria to optimize calculated with non-integer stratum sample sizes;\cr
\code{"takeall"}: the number of take-all strata for the corresponding boundaries, after making, if needed, an adjustment for non requested take-all strata;\cr 
\code{"step"}: the step in the boundary modification (for Kozak's algorithm only);\cr
\code{"iter"}: the iteration identification number;\cr 
\code{"run"}: the run identification number (for Kozak's algorithm only) because when \code{trymany=TRUE} or \code{rep}>1 the algorithm is run more than once.\cr
For Kozak's algorithm, a row is added to \code{iter.detail} only for accepted boundaries, i.e. boundaries decreasing the optimization criteria while respecting the conditions on \eqn{N_h} and \eqn{n_h} (see \bold{Details} for the algorithm description). }

  \item{niter }{Only if an iterative algorithm is run: the total number of iterations before convergence of the algorithm. For Kozak's algorithm, this number includes all the iterations, even the ones with discarded boundaries. \code{niter} is a vector if the algorithm was run more than once and if more than one run lead to the output solution. }
  
  \item{converge }{Only if an iterative algorithm is run: a logical indicating if the algorithm has converged (see \bold{Details}). }
  
  \item{run.detail }{Only for Kozak's algorithm when it is run more than once: a data frame giving the final solution for every run of the algorithm. It contains, for each run:\cr
\code{"bh"}: the stratum boundaries at convergence;\cr
\code{"opti.nh"}: the value of the criteria to optimize calculated with integer stratum sample sizes;\cr
\code{"opti.nhnonint"}: the value of the criteria to optimize calculated with non-integer stratum sample sizes;\cr
\code{"takeall"}: the number of take-all strata for the corresponding boundaries, after making, if needed, an adjustment for non requested take-all strata;\cr 
\code{"niter"}: the number of iterations of the algorithm;\cr
\code{"initbh.type"}: the type of initial boundaries requested (\code{"initbh"} for user-given boundaries, \code{"cumrootf"}, \code{"geo"} or \code{robust});\cr
\code{"ibh"}: the initial boundaries used to run the algorithm;\cr
\code{"rep"}: the repetition number (because the algorithm is run \code{rep} times for each set of initial boundaries tried).\cr
When \code{trymany=TRUE}, if a set of initial boundaries do not meet the conditions on \eqn{N_h} and \eqn{n_h} (see \bold{Details} for the algorithm description), the algorithm is not run with these boundaries and no rows are added to \code{run.detail}.}

  \item{run.min }{Only for Kozak's algorithm when it is run more than once: the identification number of every algorithm run leading to the optimal plan, i.e. the rows of \code{run.detail} containing the proposed solution. }
  
  \item{sol.detail}{Only for Kozak's algorithm when the number of possible sets of boundaries is lower than \code{minsol}: a data frame giving information for the possible solutions fulfilling the conditions \eqn{N_h \geq}{Nh >=} \code{minNh} and \code{nh} > 0. It contains:\cr
\code{"bh"}: the stratum boundaries;\cr
\code{"Nh"}: the number of units in each stratum;\cr
\code{"opti.nh"}: the value of the criteria to optimize calculated with integer stratum sample sizes;\cr
\code{"opti.nhnonint"}: the value of the criteria to optimize calculated with non-integer stratum sample sizes;\cr
\code{"takeall"}: the number of take-all strata for the corresponding boundaries, after making, if needed, an adjustment for non requested take-all strata. }

  \item{sol.min}{Only for Kozak's algorithm when the number of possible sets of boundaries is lower than \code{minsol}: the rows of \code{sol.detail} containing the optimal solution.}
  
  \item{nsol}{Only for Kozak's algorithm : the number of possible sets of boundaries (see \bold{Details}).}   
}
\references{
Baillargeon, S., Rivest, L.-P., Ferland, M. (2007). Stratification en enquetes entreprises : Une revue et quelques avancees. \emph{Proceedings of the Survey Methods Section, 2007 SSC Annual Meeting}. \url{ http://www.ssc.ca/survey/documents/SSC2007_S_Baillargeon.pdf }

Baillargeon, S. and Rivest, L.-P. (2009). A general algorithm for univariate stratification. \emph{International Stratification Review}, \bold{77}(3), 331-344.

Baillargeon, S. and  Rivest L.-P. (2011). The construction of stratified designs in R with the package stratification. \emph{Survey Methodology}, \bold{37}(1), 53-65. \url{ http://www.statcan.gc.ca/pub/12-001-x/2011001/article/11447-eng.pdf }

Dalenius, T. and Hodges, J.L., Jr. (1959). Minimum variance stratification. \emph{Journal of the American Statistical Association}, \bold{54}, 88-101.

Gunning, P. and Horgan, J.M. (2004). A new algorithm for the construction of stratum boundaries in skewed populations. \emph{Survey Methodology}, \bold{30}(2), 159-166.

Gunning, P. and Horgan, J.M. (2007). Improving the Lavallee and Hidiroglou algorithm for stratification of skewed populations. \emph{Journal of Statistical Computation and Simulation}, \bold{77}(4), 277-291.

Kozak, M. (2004). Optimal stratification using random search method in agricultural surveys. \emph{Statistics in Transition}, \bold{6}(5), 797-806.

Lavallee, P. and Hidiroglou, M.A. (1988). On the stratification of skewed populations. \emph{Survey Methodology}, \bold{14}, 33-43.

Rivest, L.-P. (1999). Stratum jumpers: Can we avoid them? \emph{Proceedings of the Section on Survey Methods Research of the American Statistical Association}, 64-72.

Rivest, L.-P. (2002). A generalization of the Lavallee and Hidiroglou algorithm for stratification in business surveys. \emph{Survey Methodology}, \bold{28}(2), 191-198.

Sethi, V. K. (1963). A note on optimum stratification of populations for estimating the population means. \emph{The Australian Journal of Statistics}, \bold{5}, 20-33.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{print.strata}}, \code{\link{plot.strata}}, \code{\link{strata.cumrootf}}, \code{\link{strata.geo}} }
\examples{
#########################################################
### Sethi's algorithm versus Kozak's algorithm

data(UScities)

# LACK OF CONVERGENCE
# Here is an example of numerical difficulties met with Sethi but not with Kozak  
Sethi <- strata.LH(x=UScities, CV=0.01, Ls=3, alloc=c(0.35,0.35,0), takenone=0, takeall=1,
         rh=1, model="loglinear", model.control=list(beta=1, sig2=0.5, ph=0.85), 
         algo="Sethi", algo.control=list(maxiter=20))
Sethi
Sethi$iter.detail[1:5,]
# Kozak's algorithm with arithmetic initial boundaries 
# (default initial boundaries for Sethi's algorithm) 
Kozak<-strata.LH(x=UScities, initbh=c(18,27), CV=0.01, Ls=3, alloc=c(0.35,0.35,0),
       takenone=0, takeall=1, rh=1, model="loglinear",
       model.control=list(beta=1, sig2=0.5, ph=0.85), algo="Kozak")
Kozak
Kozak$iter.detail[Kozak$iter.detail[,"run"]==Kozak$run.min[1],]
# Looking at the iteration history for the optimization with Sethi and Kozak, 
# we see that the initial boundaries are very close from the optimal ones. 
# Kozak reaches very quickly a minimum. However, Sethi increases n instead of 
# minimizing it and afterwards it oscillates between two sets of boundaries 
# without converging.

# LOCAL MINIMUM
# In this example, Sethi's algorithm obviously reaches a local minimum since Kozak 
# proposes a much smaller n.
Sethi<-strata.LH(x=UScities, CV=0.01, Ls=4, alloc=c(0.5,0,0), takenone=0, takeall=1, 
       rh=0.85, model="loglinear", model.control=list(beta=1.1, sig2=0, ph=1),
       algo="Sethi")
Sethi
Kozak<-strata.LH(x=UScities, CV=0.01, Ls=4, alloc=c(0.5,0,0), takenone=0, takeall=1,
       rh=0.85, model="loglinear", model.control=list(beta=1.1, sig2=0, ph=1),
       algo="Kozak")
Kozak


#########################################################
### Take-none stratum

data(MRTS)

# As illustrated in the following example (presented in Baillargeon and Rivest 2011), 
# it is sometimes beneficial to include a take-none stratum in the stratified design 
# (possibly with a bias penalty lower than 1).
notn <- strata.LH(x=MRTS, CV=0.1, Ls=3, alloc=c(0.5,0,0.5))
notn
tn1 <- strata.LH(x=MRTS, CV=0.1, Ls=3, alloc=c(0.5,0,0.5), takenone=1)
tn1
tn0.5 <- strata.LH(x=MRTS, CV=0.1, Ls=3, alloc=c(0.5,0,0.5), takenone=1, bias.penalty=0.5)
tn0.5

# Note: Sethi does not converge here. This occurs often with a take-none stratum.
tn1.Sethi <- strata.LH(x=MRTS, CV=0.1, Ls=3, alloc=c(0.5,0,0.5), takenone=1, algo="Sethi")
tn1.Sethi


#########################################################
### Automatic detection of a take-all stratum

data(UScities)

# As in the following example, a beneficial take-all stratum is not always detected 
# by the algorithm. Therefore, it is often a good idea to obtain a stratified design
# with and without a take-all stratum and to compare the results.
sans<-strata.LH(x=UScities, n=300, Ls=3, alloc=c(0.35,0.35,0), takeall=0,
      model="loglinear", model.control=list(beta=0.9, sig2=0, ph=1),
      algo.control=list(trymany=FALSE))
sans
avec<-strata.LH(x=UScities, n=300, Ls=3, alloc=c(0.35,0.35,0), takeall=1,
      model="loglinear", model.control=list(beta=0.9,sig2=0,ph=1),
      algo.control=list(trymany=FALSE))
avec
# We see that for the target sample size, the anticipated CV is 17\% lower with a 
# take-all stratum (0.01081053 vs 0.009002313).


#########################################################
### Models for the discrepancy between Y and X

data(Sweden)

# LOGLINEAR MODEL WITH MORTALITY: see help(Sweden)

# HETEROSCEDASTIC LINEAR MODEL:  We fix gamma=2.
beta.lin <- mean(Sweden$RMT85/Sweden$REV84)
sig2.lin <- var(Sweden$RMT85/Sweden$REV84)
strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
          model="linear", model.control=list(beta=beta.lin, sig2=sig2.lin, gamma=2),
          algo="Kozak")
# Verification of equation 3.6 of Rivest (2002)
beta.log <- 1
sig2.log <- log(1+sig2.lin/beta.lin^2)
strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
          model="loglinear", model.control=list(beta=beta.log, sig2=sig2.log, ph=1),
          algo="Kozak")
# The two models give the same stratified design.

# RANDOM REPLACEMENT MODEL: example in Rivest (1999)
strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
       model="none", algo="Sethi")     # Table 1 with a different rounding of the nh's
e0 <- strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
       model="none", algo="Kozak")
e0                                    # Better than Sethi
var.strata(e0, y=Sweden$RMT85)
e0.001 <- strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
          model="random", model.control=list(epsilon=0.011), algo="Kozak")
e0.001                          # Table 2 part 3 with a different rounding of the nh's
var.strata(e0.001 ,y=Sweden$RMT85)

}
\keyword{ survey }
