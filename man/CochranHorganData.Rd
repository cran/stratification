\name{CochranHorganData}
\alias{Debtors}
\alias{UScities}
\alias{UScolleges}
\alias{USbanks}
\docType{data}
\title{ Populations Analyzed in Gunning and Horgan (2004) and Cochran (1961) }
\description{
The first population \code{Debtors} is an accounting population of debtors in an Irish firm, detailed in Horgan (2003). The other populations are three of the skewed populations in Cochran (1961). These are:\cr
\code{UScities}: the population in thousands of US cities in 1940;\cr
\code{UScolleges}: the number of students in four-year US colleges in 1952-1953;\cr
\code{USbanks}: the resources in millions of dollars of large commercial US banks.
}
\usage{
Debtors
UScities
UScolleges
USbanks
}
\format{
  The formats of these data sets are, respectively:\cr
 num [1:3369] 40 40 40 40 40 40 40 40 40 40 ...\cr
 num [1:1038] 10 10 10 10 10 10 10 10 10 10 ...\cr
 num [1:677] 200 201 202 202 207 210 211 213 215 217 ...\cr
 num [1:357] 70 71 72 72 72 73 73 73 73 73 ...\cr
}
\source{
Jane M. Horgan
}
\references{
Cochran, W.G. (1961). Comparison of methods for determining stratum boundaries. \emph{Bulletin of the International Statistical Institute}, \bold{32}(2), 345-358.

Gunning, P. and Horgan, J.M. (2004). A new algorithm for the construction of stratum boundaries in skewed populations. \emph{Survey Methodology}, \bold{30}(2), 159-166.

Horgan, J.M. (2003). A list sequential sampling scheme with applications in financial auditing. \emph{IMA Journal of Management Mathematics}, \bold{14}, 1-18.
}
\examples{
### Reproduction of the results in Table 4 and Table 7 part 3 (case L=5) of  
### Gunning and Horgan (2004). The differences in the nh come from different  
### rounding. The more important differences observed for the cumulative  
### root frequency method are due to the use of different numbers of classes.
strata.geo(x=Debtors, n=100, Ls=5, alloc=c(0.5,0,0.5))
strata.cumrootf(x=Debtors, n=100, Ls=5, alloc=c(0.5,0,0.5), nclass=40)
strata.LH(x=Debtors, CV=0.0360, Ls=5, alloc=c(0.5,0,0.5), takeall=1, algo="Sethi")

strata.geo(x=UScities, n=100, Ls=5, alloc=c(0.5,0,0.5))
strata.cumrootf(x=UScities, n=100, Ls=5, alloc=c(0.5,0,0.5), nclass=40)
strata.LH(x=UScities, CV=0.0144, Ls=5, alloc=c(0.5,0,0.5), takeall=1, algo="Sethi")

strata.geo(x=UScolleges, n=100, Ls=5, alloc=c(0.5,0,0.5))
strata.cumrootf(x=UScolleges, n=100, Ls=5, alloc=c(0.5,0,0.5), nclass=40)
strata.LH(x=UScolleges, CV=0.0184, Ls=5, alloc=c(0.5,0,0.5), takeall=1, algo="Sethi")

strata.geo(x=USbanks, n=100, Ls=5, alloc=c(0.5,0,0.5))
strata.cumrootf(x=USbanks, n=100, Ls=5, alloc=c(0.5,0,0.5), nclass=40)
strata.LH(x=USbanks, CV=0.0110, Ls=5, alloc=c(0.5,0,0.5), takeall=1, algo="Sethi")
}
\keyword{datasets}
