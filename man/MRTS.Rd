\name{MRTS}
\alias{MRTS}
\docType{data}
\title{ Simulated Data from the Monthly Retail Trade Survey (MRTS) of Statistics Canada }
\description{
This data set is a vector containing the simulated values of a realistic stratification variable: the size measure used for Canadian retailers in the Monthly Retail Trade Survey (MRTS) carried out by Statistics Canada. This size measure is created using a combination of independent survey data and three administrative variables from the corporation tax return. The MRTS aims at estimating sales from retailers, which are a key monthly indicator of consumer purchasing patterns in Canada.  
}
\usage{data(MRTS)}
\format{
  The format is:
 num [1:2000] 141 209 238 257 261 ...
}
\source{
The data set has been simulated with the command:\cr
\code{exp(rst(n=2000, location=9.7093354, scale=0.7885551,}\cr
\code{shape=-0.6384867, df=5.6243544))}.\cr
The function \code{rst} comes from the package \pkg{sn}. It generates random numbers for the skew-t distribution. The parameters of the distribution have been estimated with the assistance of Michel Ferland from Statistics Canada to be representative of the measure of size used in the MRTS.
}
\references{
Baillargeon, S., Rivest, L.-P., Ferland, M. (2007). Stratification en enquetes entreprises : Une revue et quelques avancees. \emph{Proceedings of the Survey Methods Section, 2007 SSC Annual Meeting}. \url{ http://www.ssc.ca/survey/documents/SSC2007_S_Baillargeon.pdf }

Information on Statistics Canada Monthly Retail Trade Survey: \url{http://www.statcan.ca/cgi-bin/imdb/p2SV.pl?Function=getSurvey&SDDS=2406&lang=en&db=IMDB&dbg=f&adm=8&dis=2}
}
\examples{
# Production of results similar to those in Table 1 of Baillargeon, Rivest
# and Ferland (2007). The results are not the same because calculations in
# the paper were conducted on real data whereas, for confidentiality reason, 
# the MRTS data included in the package is simulated.
data(MRTS)
geo<-strata.geo(x=MRTS,CV=0.01,Ls=4,alloc=c(0.5,0,0.5))
geo
aRRMSE.geo<-var.strata(geo,model="loglinear",
            model.control=list(beta=0.9,sig2=0.015,ph=c(0.8,0.9,0.95,1)))
aRRMSE.geo$RRMSE
plot(geo,logscale=TRUE)
# The geometric method does not perform well because of some small units

cumrootf<-strata.cumrootf(x=MRTS,nclass=500,CV=0.01,Ls=4,alloc=c(0.5,0,0.5))
cumrootf
aRRMSE.cum<-var.strata(cumrootf,,rh=c(0.85,0.9,0.9,1),model="loglinear",
            model.control=list(beta=0.9,sig2=0.015,ph=c(0.8,0.9,0.95,1)))
aRRMSE.cum$RRMSE

LH<-strata.LH(x=MRTS,CV=0.01,Ls=4,alloc=c(0.5,0,0.5),takeall=1,algo="Sethi")
LH
aRRMSE.LH<-var.strata(LH,rh=c(0.85,0.9,0.9,1),model="loglinear",
           model.control=list(beta=0.9,sig2=0.015,ph=c(0.8,0.9,0.95,1)))
aRRMSE.LH$RRMSE

LH.full<-strata.LH(x=MRTS,CV=0.01,Ls=4,alloc=c(0.5,0,0.5),takeall=1,
           algo="Sethi",rh=c(0.85,0.9,0.9,1),model="loglinear",
           model.control=list(beta=0.9,sig2=0.015,ph=c(0.8,0.9,0.95,1)))
LH.full
aRRMSE.LH.full<-var.strata(LH.full,rh=c(0.85,0.9,0.9,1),model="loglinear",
           model.control=list(beta=0.9,sig2=0.015,ph=c(0.8,0.9,0.95,1)))
aRRMSE.LH.full$RRMSE
}
\keyword{datasets}
