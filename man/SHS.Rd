\name{SHS}
\alias{SHS}
\docType{data}
\title{ Data from the 2001 Survey of Household Spending (SHS), Statistics Canada }
\description{
This data set contains some variables from the 2001 Survey of Household Spending (SHS) carried out by Statistics Canada. The main purpose of this survey is to obtain detailed information about household spending during the reference year.
}
\usage{data(SHS)}
\format{
  A data frame with 16057 observations on the following 7 variables.
  \describe{
    \item{\code{CASEID}}{Identification number}
    \item{\code{WEIGHT}}{Weight at household level}
    \item{\code{PROVINCP}}{Province or territory code}
    \item{\code{URBRUR}}{Urban rural code}
    \item{\code{URBSIZEP}}{Size of area of residence code}
    \item{\code{HHINCTOT}}{Household income before taxes}
    \item{\code{M101}}{Household spending on recreation}
  }
}
\details{
In this package, \code{HHINCTOT} is used as a stratification variable and \code{M101} as a survey variable.
}
\source{
Income Statistics Division, Statistics Canada.
}
\references{
Information on Statistics Canada Survey of Household Spending: \url{http://www.statcan.ca/cgi-bin/imdb/p2SV.pl?Function=getSurvey&SDDS=3508&lang=en&db=IMDB&dbg=f&adm=8&dis=2}
}
\examples{
data(SHS)

# Estimation of the response rate
X<-SHS$HHINCTOT[SHS$HHINCTOT>0]
Y<-SHS$M101[SHS$HHINCTOT>0]
Y[Y<0]<-0
p<-sum(Y>0)/length(Y)

# Study of the relationship between X and Y for the active units
Xactive<-SHS$HHINCTOT[SHS$HHINCTOT>0&SHS$M101>0]
Yactive<-SHS$M101[SHS$HHINCTOT>0&SHS$M101>0]
plot(log(Xactive),log(Yactive))
# Extreme values are omitted for a more robust estimation
keep<-Xactive/Yactive>quantile(Xactive/Yactive,0.01)&
      Xactive/Yactive<quantile(Xactive/Yactive,0.99)    
plot(log(Xactive)[keep],log(Yactive)[keep])
reg<-lm( log(Yactive)[keep]~log(Xactive)[keep] )
summary(reg)

# Stratification assuming X=Y
nomodel<-strata.LH(x=X,CV=0.05,Ls=3,alloc=c(0.5,0,0.5),takeall=0,model="none")
nomodel
var.strata(nomodel,y=Y) # The target CV is far from being reached

# Stratification taking into account a loglinear model with mortality 
# between X and Y, using the estimated parameters values
model<-strata.LH(x=X,CV=0.05,Ls=3,alloc=c(0.5,0,0.5),takeall=0,model="loglinear",
       model.control=list(beta=reg$coef[2],sig2=summary(reg)$sigma^2,ph=0.97),
       algo.control=list(rep=2))
model
var.strata(model,y=Y) # The target CV is reached
}
\keyword{datasets}
