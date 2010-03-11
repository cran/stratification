\name{Sweden}
\alias{Sweden}
\docType{data}
\title{ The MU284 Population of Sweden Municipalities from Sarndal et al. (1992) }
\description{
This data set comes from Sarndal et al.'s book (1992), Appendix B. It contains different variables that describe 284 municipalities in Sweden.
}
\usage{data(Sweden)}
\format{
  A data frame with 284 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{Identifier running from 1 to 284}
    \item{\code{P85}}{1985 population (in thousands)}
    \item{\code{P75}}{1975 population (in thousands)}
    \item{\code{RMT85}}{Revenues from the 1985 municipal taxation (in millions of kronor)}
    \item{\code{CS82}}{Number of Conservative seats in municipal council}
    \item{\code{SS82}}{Number of Social-Democratic seats in municipal council}
    \item{\code{S82}}{Total number of seats in municipal council}
    \item{\code{ME84}}{Number of municipal employees in 1984}
    \item{\code{REV84}}{Real estate values according to 1984 assessment (in millions of kronor)}
    \item{\code{REG}}{Geographic region indicator}
    \item{\code{CL}}{Cluster indicator (a cluster consists of a set of neighboring municipalities)}
  }
}
\details{
In this package, \code{REV84} is used as a stratification variable and \code{RMT85} as a survey variable.
}
\source{
Sarndal, C. E., Swensson, B. and Wretman, J. (1992). \emph{Model Assisted Survey Sampling}. Springer Verlag, New York.
}
\references{
Rivest, L.-P. (2002). A generalization of the Lavallee and Hidiroglou algorithm for stratification in business surveys. \emph{Survey Methodology}, \bold{28}(2), 191-198.
}
\examples{
data(Sweden)
X<-Sweden$REV84
Y<-Sweden$RMT85

# Study of the relationship between X and Y
plot(log(X),log(Y))
# Extreme values are omitted for a more robust estimation
keep<-X/Y>quantile(X/Y,0.03)&X/Y<quantile(X/Y,0.97)
plot(log(X)[keep],log(Y)[keep])
reg<-lm( log(Y)[keep]~log(X)[keep] )
summary(reg)

# Stratification assuming X=Y
nomodel<-strata.LH(x=X,CV=0.05,Ls=3,alloc=c(0.5,0,0.5),takeall=1,model="none")
nomodel
var.strata(nomodel,y=Y) # The target CV is not reached

# Stratification taking into account a loglinear model between X and Y, 
# using the estimated parameters values
model<-strata.LH(x=X,CV=0.05,Ls=3,alloc=c(0.5,0,0.5),takeall=1,model="loglinear",
       model.control=list(beta=reg$coef[2],sig2=summary(reg)$sigma^2,ph=1))
model
var.strata(model,y=Y) # The target CV is reached
}
\keyword{datasets}
