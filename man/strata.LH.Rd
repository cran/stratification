\name{strata.LH}
\alias{strata.LH}
\title{ Generalized Lavallee-Hidiroglou Method of Strata Construction }
\description{
This function constructs optimal strata with a generalized Lavallee-Hidiroglou (1998) method.
The function uses Kozak's (2004) algorithm by default, but it can also apply Sethi's (1963) algorithm (argument \code{algo="Sethi"}).
The function takes into account a discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can consider a loglinear model with mortality between the variables. With Kozak's algorithm, two additional models are implemented: an heteroscedastic linear model and a random replacement model as in Rivest (2002). The determination of the optimal boundaries also incorporates, if desired, an anticipated non-response, a take-all stratum for the large units and a take-none stratum for the small units.\cr\cr
Sethi's algorithm is not used by default because it can be numerically unstable, especially with a take-none stratum. Better results were obtained with Kozak's algorithm in our numerical experiments.
}
\usage{
strata.LH(x, initbh, n = NULL, CV = NULL, Ls = 3, certain = NULL,
          alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), takenone = 0, 
          bias.penalty = 1, takeall = 0, rh = rep(1, Ls), 
          model = c("none", "loglinear", "linear", "random"), 
          model.control = list(), algo = c("Kozak", "Sethi"), 
          algo.control = list())
}
\arguments{
  \item{x}{ A vector containing the values of the stratification variable \eqn{X} for every unit in the population. }
  \item{initbh}{ A vector of \code{Ls}-1, or \code{Ls} if \code{takenone}=1, initial stratum boundaries  (see \bold{Details}).}  
  \item{n}{ A numeric: the target sample size. It has no default value. The argument \code{n} or the argument \code{CV} must be input. }
  \item{CV}{ A numeric: the target coefficient of variation or relative root mean squared error if \code{takenone}=1. It has no default value. The argument \code{CV} or the argument \code{n} must be input. }
  \item{Ls}{ A numeric: the number of sampled strata (excluding the take-none stratum). The default is 3. }
  \item{certain}{ A vector giving the position, in the vector \code{x}, of the units that must be included in the sample (see \code{\link{stratification-package}}). By default \code{certain} is \code{NULL}, which means that no units are chosen a priori to be in the sample. }
  \item{alloc}{ A list specifying the allocation scheme. The list must contain 3 numerics for the 3 exponents \code{q1}, \code{q2} and \code{q3} in the general allocation scheme (see \code{\link{stratification-package}}). The default is Neyman allocation (\code{q1}=\code{q3}=0.5 and \code{q2}=0)}
  \item{takenone}{ A numeric: the number of take-none strata (0 or 1). The default is 0, i.e. no take-none stratum is included.}
  \item{bias.penalty}{ A numeric between 0 and 1 giving the penalty for the bias in the anticipated mean squared error (MSE) of the survey estimator (see \code{\link{stratification-package}}). This argument is relevant only if \code{takenone}=1. The default is 1.}
  \item{takeall}{ A numeric: the number of take-all strata. The default is 0, i.e. no take-all stratum is included.}
  \item{rh}{ A vector giving the anticipated response rates in each of the \code{Ls} sampled strata. A single number can be given if the rates do not vary among strata. The default is 1 in each stratum. }
  \item{model}{ A character string identifying the model used to describe the discrepancy between the stratification variable \eqn{X} and the survey variable \eqn{Y}. It can be \code{"none"} if one assumes \eqn{Y=X}, \code{"loglinear"} for the loglinear model with mortality, \code{"linear"} for the heteroscedastic linear model or \code{"random"} for the random replacement model (see \code{\link{stratification-package}} for a description of these models). The default is \code{"none"}, so the original Lavallee-Hidiroglou (1998) method of strata construction is used. The last two models \code{"linear"} and \code{"random"} are only available with Kozak's algorithm.}
  \item{model.control}{ A list of model parameters (see \code{\link{stratification-package}}). The default values of the parameters correspond to the model \eqn{Y=X}. }
  \item{algo}{ A character string identifying which optimization algorithm is to be used. It can be "Kozak" for Kozak's (2004) algorithm or "Sethi" for Sethi's (1963) algorithm. The default is "Kozak". }
  \item{algo.control}{ A list of parameters to control the optimization algorithm (see \bold{Details}).}
}
\details{
INITIAL STRATUM BOUNDARIES\cr
%The default value of \code{initbh} is the set of \code{Ls}-1 equidistant points along the range of the \eqn{X}-values 
%(arithmetic starting points of Gunning and Horgan (2007)). Let's define \eqn{d=(max(X)-min(X))/Ls}, the default initial
% boundaries are \eqn{b_h = min(X) + h \times d}{bh = min(X) + h x d} for \eqn{h=1,\ldots,Ls-1}. 
The default initial boundaries \code{initbh} are obtained with the function \code{\link{strata.cumrootf}}
(cumulative root frequency method). 

If \code{takenone}=1 and \code{initbh} is of size \code{Ls}-1, the initial boundary of the take-none stratum is set to 
the first percentile of \code{x}. If this first percentile is equal to the minimum value of \code{x}, this initial 
boundary would lead to an empty take-none stratum. In that case, the initial boundary of the take-none stratum is 
rather set to the second smallest value of \code{x}.

Sometimes, the specified initial boundaries are not suitable for Kozak's algorithm. As written below, this algorithm 
verifies at each iteration that the sampled strata contains at least \code{minNh} units and that they have positive 
sample sizes \eqn{n_h}{nh}. If the initial boundaries do not meet these conditions, the algorithm might not be able 
to change the boundaries at all because any change would be rejected for failing to meet the conditions. Therefore, 
in such a situation, the specified initial boundaries are replaced by robust boundaries. The robust boundaries give 
an empty take-none stratum if such a stratum is requested, take-all strata as small as possible, and take-some strata 
with approximately the same number of unique \eqn{X}-values. The sampled strata given by the robust boundaries contain 
at least \code{minNh} units. These robust boundaries give, in most situation, a non-negative optimization criteria 
\eqn{n} if a target \code{CV} is requested (see Rivest and Baillargeon 2009 for a discussion of possible negative 
\eqn{n} values).\cr  

ALGORITHMS\cr
\bold{Sethi:} The formulas implemented for Sethi's algorithm are presented in Baillargeon, Rivest and Ferland (2007).\cr
\bold{Original Kozak:} The original Kozak's algorithm is described in Kozak (2004). At each iteration, it chooses a stratum boundary at random and a random modification for this boundary among the 2*\code{maxstep} possible alternatives. If this modification reduces the optimization criteria, creates sampled strata containing at least \code{minNh} units and leads to positive \eqn{n_h}{nh}, it is accepted. Otherwise, the boundaries are not changed. If the boundaries remain unchanged for \code{maxstill} consecutive iterations, the algorithm stops. This algorithm is run \code{rep} times. Note: Kozak's (2004) used the condition \eqn{n_h \geq 2}{nh >= 2} instead of \eqn{n_h > 0}{nh > 0}. We chose to modify this requirement because we noticed that more severe conditions sometimes prevented the algorithm from selecting a path leading to the optimal solution. According to our numerical experiments, if one sets a posteriori to 2 the \eqn{n_h}{nh}'s that are equal to 1 at the end of the algorithm, the new sample size is smaller than or equal to the one obtained with the condition \eqn{n_h \geq 2}{nh >= 2}.\cr
\bold{Modified Kozak:} The modified Kozak's algorithm is a non-random version of Kozak's algorithm. At each iteration of the algorithm, the (\eqn{L-1})*2*\code{maxstep} possible boundaries modifications among which the original Kozak's algorithm chooses at random are tried and the modification giving the largest decrease in the optimization criteria is kept. The modification must also create sampled strata containing at least \code{minNh} units and it must lead to positive \eqn{n_h}{nh}. This algorithm is slower than the original Kozak's algorithm.\cr

For data sets with a small number of unique \eqn{X}-values \code{Nu}, Kozak's algorithm is often caught in a local minimum. In such situations, a complete enumeration of all possible sets of boundaries is feasible. When \code{takenone=0}, there are \code{choose(Nu-1,L-1)} possible sets of boundaries. The function \code{strata.LH} carries out a complete enumeration of all possible cases if the number of possible sets of boundaries is lower than the \code{algo.control} argument \code{minsol}.\cr 

ALGORITHM PARAMETERS\cr
The \code{algo.control} argument is a list to supply any of the following parameters which control the algorithm. Sethi's algorithm only uses the first argument \code{maxiter}.
\describe{
  \item{\code{maxiter}}{ A numeric: the maximal number of iterations. The default is 500 for Sethi's algorithm, 10 000 for the 
                         original Kozak's algorithm and 3000 for the modified Kozak's algorithm.  }
  \item{\code{method}}{ A character string identifying which Kozak's algorithm is to be used. It can be "original" for the original 
                        Kozak's algorithm or "modified" for the modified Kozak's algorithm. The default is "original".}
  \item{\code{minsol}}{ A numeric for Kozak's algorithms only: the minimum number of possible solutions below which the algorithm 
                        is not run. If the number of possible sets of boundaries is lower than \code{minsol}, every set of boundaries 
                        is tried. This ensures that the global minimum is reached. The default value of \code{minsol} is 1000. This 
                        parameter has to take a value between 2 and 1 000 000.\cr }
  \item{\code{idopti}}{ A character string for Kozak's algorithms only: this argument determines which stratum sample sizes are going 
                        to be used to calculate the optimization criteria. It can take the value \code{"nh"} (criteria calculated with 
                        integer, so rounded, \eqn{n_h}{nh}) or \code{"nhnonint"} (criteria calculated with non integer \eqn{n_h}{nh}).
                        The default value is \code{"nh"}. When a complete enumeration is performed, the criteria is calculated 
                        with integer stratum sample sizes since they are the one used in practice. For Sethi's algorithm, the 
                        criteria can only be calculated with non integer stratum sample sizes. \cr }
  \item{\code{minNh}}{ A numeric for Kozak's algorithms only: the minimum number of units required in each sampled stratum (no 
                       restriction is put on the take-none stratum, if included). \code{minNh} must be greater or equal to 2, 
                       which is the default. }
  \item{\code{maxstep}}{ A numeric for Kozak's algorithms only: the maximal step. The default is 3 for the modified Kozak's algorithm
                         and \code{pmin(ceiling(N1/15),75)} for the original Kozak's algorithm, where N1 stands for the number of 
                         unique values in the input data. }
  \item{\code{maxstill}}{ A numeric for the original Kozak's algorithm only: the maximal number of iterations without a change in 
                          the boundaries. The default is \code{floor(maxstep*100/3)}.}
  \item{\code{rep}}{ A numeric for Kozak's original algorithm only: the number of repetitions of the algorithm.
                     The default is 5. 
                     This argument can be given the special value \code{"change"}. 
                     If \code{rep="change"} the algorithm is repeated 
                     3 times for 9 combinations of initial boundaries and \code{maxstep} values. The initial boundaries used are the 
                     cumulative root frequency boundaries, the geometric boundaries or the default initial boundaries. The 
                     \code{maxstep} values used are 3, \eqn{min( 10,}number of unique \eqn{X}-values) or 
                     \eqn{min(}\code{ceiling}(number of unique \eqn{X}-values/10),100). The \code{maxstill} argument is given a 
                     value proportional to \code{maxstep}. The option \code{rep="change"} can be long to run, but it increases the 
                     chances for the algorithm to reach the global minimum.}
}

CONVERGENCE\cr
It is possible for the algorithm not to converge. In this case, a warning is printed. The only possible cause of non-convergence for Kozak's algorithm, original or modified, is to reach the maximum number of iterations before the stopping rule has been met. For the algorithm to converge, the argument \code{maxiter} has to be increased. On the other hand, non-convergence of Sethi's algorithm has several possible causes:
\itemize{
\item a division by zero caused by an empty stratum can occur ;
\item a division by zero caused by a 0 stratum variance can occur ;
\item the square root of a negative number (negative discriminant) can occur ;
\item the maximum number of iterations can be reached (often because the algorithm is caught in a loop of non-optimal sets of boundaries).
}
If a non-convergence happens, the user can try to change the initial boundaries or the model parameters. The user can also choose to work with Kozak's algorithm which should converge given an appropriate maximum number of iterations.\cr

Let's note however that even if the algorithm converges, it is not guaranteed that it has reached a global minimum. Local minimum can occur with both Sethi's and Kozak's algorithms (Rivest and Baillargeon, 2009).
}
\value{
  \item{Nh }{A vector of length \eqn{L} containing the population sizes \eqn{N_h}{Nh}, i.e. the number of units in each stratum.}
  \item{nh }{A vector of length \eqn{L} containing the sample sizes \eqn{n_h}{nh}, i.e. the number of units to sample in each stratum. See \code{\link{stratification-package}} for information about the rounding used to get these values.}
  \item{n }{The total sample size (\code{sum(nh)}).}
  \item{nh.nonint }{A vector of length \eqn{L} containing the non-integer values of the sample sizes \eqn{n_h}{nh}.}
  \item{certain.info }{A vector giving statistics for the certain stratum (see \code{\link{stratification-package}}). It contains \code{Nc}, the number of units chosen a priori to be in the sample, and \code{meanc}, the anticipated mean of \eqn{Y} for these units. }
  \item{opti.criteria }{The final value of the criteria to optimize : either the sum of the stratum sample sizes (if a target 
                        \code{CV} was given) or the RRMSE associated to the stratum sample sizes (if a target \code{n} was given). 
                        If the algorithm parameter \code{idopti} takes the value \code{"nh"}, the stratum sample sizes used for the 
                        computation of this criteria are the integer (rounded) ones. The non-integer stratum sample sizes are used 
                        if \code{idopti="nhnonint"} (the default).}
  \item{bh }{A vector of the \eqn{L-1} optimal stratum boundaries found by the algorithm.}
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
  \item{initbh }{Only if an iterative algorithm is run: a vector of the \eqn{L-1} initial stratum boundaries. Most of the time, this will be a copy of the input initial boundaries, but sometimes these input boundaries have to be replaced by robust boundaries (see \bold{Details}). }
  \item{iter.detail }{Only if an iterative algorithm is run. A matrix giving the detail of every iteration of the algorithm. Because of the automatic adjustment for a take-all stratum, the algorithm can be run more than once by adding 1 to the \code{takeall} argument each time. All the runs of the algorithm are included in \code{iter.detail}. The matrix \code{iter.detail} contains, for each iteration, the stratum boundaries, the step for Kozak's algorithm, the value of the criteria to optimize calculated with the integer and non integer stratum sample sizes, the iteration identification number and the run identification number. For the original Kozak's algorithm, the repetition number is put in \code{iter.detail} last column. Also, for this algorithm, a row is added to \code{iter.detail} only for iterations that decrease the optimization criteria. }
  \item{niter }{Only if an iterative algorithm is run: the total number of iterations in the algorithm's run leading to the output solution. This number includes all the iterations, even the discarded ones. }
  \item{converge }{Only if an iterative algorithm is run: a logical indicating if the algorithm has converged (see \bold{Details}). }
  \item{rep.detail }{Only if \code{rep} is not equal to 1 with Kozak's original algorithm: a matrix giving the final solution for every repetition of the algorithm. It contains, for each repetition, the value of the criteria to optimize calculated with the integer (\code{"opti.nh"}) and non integer stratum sample sizes (\code{"opti.nhnonint"}), the number of take-all strata (\code{"takeall"}), the stratum boundaries (\code{"bh"}), the number of iterations of the algorithm for the final number of take-all strata (\code{"niter"}), the maxstep value (\code{"maxstep"}), the maxstill value (\code{"maxstill"}), the type of initial boundaries requested (\code{"initbh.type"}), the initial boundaries requested (\code{"initbh"}), a logical indicating whether the initial boundaries meet the condition \eqn{N_h \geq}{Nh >=} \code{minNh} (\code{"Nhok.initbh"}), a logical indicating whether the initial boundaries meet the condition \eqn{n_h > 0}{nh >= 0} (\code{"nhok.initbh"}), the initial boundaries used to run the algorithm (\code{"bhi"}). The boundaries \code{"bhi"} are the same as the boundaries \code{"initbh"} if \code{"Nhok.initbh"} and \code{"nhok.initbh"} are TRUE. Otherwise, \code{"bhi"} are robust boundaries (see \bold{Details}). If robust boundaries have to be used, the last two columns of \code{rep.detail} indicate whether these boundaries meet the condition \eqn{N_h \geq}{Nh >=} \code{minNh} (\code{"Nhok.robust"}) and the condition \eqn{n_h > 0}{nh >= 0} (\code{"nhok.robust"}). }
  \item{rep.min }{Only if \code{rep} is not equal to 1 with Kozak's original algorithm: the identification number of the algorithm repetition leading to the optimal plan, i.e. the row of \code{rep.detail} containing the proposed solution. }
  \item{sol.detail}{Only if the number of possible sets of boundaries is lower than \code{minsol}: a matrix giving information for the possible solutions fulfilling the conditions \eqn{N_h \geq}{Nh >=} \code{minNh} and \eqn{n_h > 0}{nh >= 0}. It contains the boundaries (\code{"bh"}), the number of units in each stratum (\code{"Nh"}), the stratum sample size (\code{"nh"}), the value of the criteria to optimize calculated with the integer (\code{"opti.nh"}) and non integer stratum sample sizes (\code{"opti.nhnonint"}), and the number of take-all strata (\code{"takeall"}). }
\item{sol.min}{Only if the number of possible sets of boundaries is lower than \code{minsol}: the row of \code{sol.detail} containing the optimal solution.}
\item{nsol}{Only if the number of possible sets of boundaries is lower than \code{minsol}: the number of possible sets of boundaries (\code{choose(Nu-1,L-1)} if \code{takenone=0}, \code{choose(Nu,L-1)} otherwise).}   
}
\references{
Baillargeon, S., Rivest, L.-P., Ferland, M. (2007). Stratification en enquetes entreprises : Une revue et quelques avancees. \emph{Proceedings of the Survey Methods Section, 2007 SSC Annual Meeting}. \url{ http://www.ssc.ca/survey/documents/SSC2007_S_Baillargeon.pdf }

Baillargeon, S. and Rivest, L.-P. (2009). A general algorithm for univariate stratification. \emph{International Stratification Review}, \bold{77}(3), 331-344.

Baillargeon, S. and  Rivest L.-P. (2011). The construction of stratified designs in R with the package stratification. \emph{Survey Methodology}, \bold{37}(1), 53-65. \url{ http://www.statcan.gc.ca/pub/12-001-x/2011001/article/11447-eng.pdf }

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

# LACK OF CONVERGENCE
# Here is an example of numerical difficulties met with Sethi but not with Kozak  
data(UScities)
Sethi<-strata.LH(x=UScities, CV=0.01, Ls=3, alloc=c(0.35,0.35,0), takenone=0, takeall=1,
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
Kozak$iter.detail[Kozak$iter.detail[,"rep"]==Kozak$rep.min,]
# Looking at the iteration history for the optimization with Sethi and Kozak, 
# we see that the initial boundaries are very close from the optimal ones. 
# Kozak reaches very quickly a minimum. However, Sethi increases n instead of 
# minimizing it and afterwards it oscillates between two sets of boundaries 
# without converging.

# LOCAL MINIMUM
# In this example, Sethi's algorithm obviously reaches a local minimum since Kozak 
# proposes a much smaller n.
data(UScities)
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

# As illustrated in the following example, it is sometimes beneficial to include a
# take-none stratum in the stratified design (possibly with a bias penalty lower than 1).
data(Debtors)
sans<-strata.LH(x=Debtors, CV=0.05, Ls=3, alloc=c(0.35,0.35,0), takenone=0,
      takeall=1, rh=c(0.8,0.9,1), model="loglinear", 
      model.control=list(beta=1.1, sig2=0, ph=1, ptakenone=1), algo="Kozak",
      algo.control=list(maxstep=50, maxstill=1000, rep=3))
sans
avec<-strata.LH(x=Debtors, initbh=sans$bh, CV=0.05, Ls=3, alloc=c(0.35,0.35,0), 
      takenone=1, takeall=1, rh=c(0.8,0.9,1), model="loglinear",
      model.control=list(beta=1.1, sig2=0, ph=1, ptakenone=1), algo="Kozak",
      algo.control=list(maxstep=25, maxstill=500, rep=2))
avec
avec0.5<-strata.LH(x=Debtors, initbh=avec$bh, CV=0.05, Ls=3, alloc=c(0.35,0.35,0),
      takenone=1, bias.penalty=0.5, takeall=1, rh=c(0.8,0.9,1), model="loglinear",
      model.control=list(beta=1.1, sig2=0, ph=1, ptakenone=1), algo="Kozak",
      algo.control=list(maxstep=25, maxstill=500, rep=2))
avec0.5

# Note: Sethi does not converge here. This occurs often with a take-none stratum.
avec.Sethi<-strata.LH(x=Debtors,CV=0.05,Ls=3,alloc=c(0.35,0.35,0),takenone=1,
      takeall=1,rh=c(0.8,0.9,1),model="loglinear",model.control=list(beta=1.1,
      sig2=0,ph=1,ptakenone=1),algo="Sethi")
avec.Sethi


#########################################################
### Automatic detection of a take-all stratum

# As in the following example, a beneficial take-all stratum is not always 
# detected by the algorithm. Therefore, it is often a good idea to obtain a
# stratified design with and without a take-all stratum and to compare 
# the results.
data(UScities)
sans<-strata.LH(x=UScities, CV=0.01, Ls=3, alloc=c(0.35,0.35,0), takenone=0, takeall=0,
      rh=1, model="loglinear", model.control=list(beta=0.9, sig2=0, ph=1), algo="Kozak")
sans
avec<-strata.LH(x=UScities, CV=0.01, Ls=3, alloc=c(0.35,0.35,0), takenone=0, takeall=1,
      rh=1, model="loglinear", model.control=list(beta=0.9,sig2=0,ph=1), algo="Kozak")
avec


#########################################################
### Models for the discrepancy between Y and X

data(Sweden)

# LOGLINEAR MODEL WITH MORTALITY: see help(Sweden)

# HETEROSCEDASTIC LINEAR MODEL:  We fix gamma=2.
beta.lin <- mean(Sweden$RMT85/Sweden$REV84)
sig2.lin <- var(Sweden$RMT85/Sweden$REV84)
strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
       model="linear", model.control=list(beta=beta.lin, sig2=sig2.lin, gamma=2),
       algo="Kozak", algo.control=list(method="modified"))
# Verification of equation 3.6 of Rivest (2002)
beta.log <- 1
sig2.log <- log(1+sig2.lin/beta.lin^2)
strata.LH(x=Sweden$REV84, CV=0.05, Ls=5, alloc=c(0.5,0,0.5), takeall=1,
       model="loglinear", model.control=list(beta=beta.log, sig2=sig2.log, ph=1),
       algo="Kozak", algo.control=list(method="modified"))
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
e0.001                         # Table 2 part 3 with a different rounding of the nh's
var.strata(e0.001 ,y=Sweden$RMT85)

}
\keyword{ survey }
