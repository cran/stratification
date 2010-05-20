\name{strata.tool}
\alias{print.strata}
\alias{plot.strata}
\title{ Functions to Visualize Stratified Designs }
\description{
\code{print.strata} prints a "strata" object, presenting the stratification information into a table.\cr

\code{plot.strata} produces a histogram of the stratification variable \eqn{X} in which the stratification boundaries are drawn. A table with the \code{Nh} and \code{nh} values is also added at the top of the plot.\cr
}
\usage{
\method{print}{strata}(x, \dots)

\method{plot}{strata}(x, logscale = FALSE, drop = 0, main = 
   paste("Graphical Representation of the Stratified Design", xname), 
   xlab, \dots)
}
\arguments{
  \item{x}{ An object of class "strata" to print or to plot. }
  \item{logscale}{ A logical indicating whether the \eqn{X} axis should be represented on the log scale or not. The default is FALSE.}
  \item{drop}{ A integer indicating how many of the largest values of \eqn{X} should be omitted in the plot. This argument is useful when some large values of \eqn{X} stretch the \eqn{X} range too much. }
  \item{main}{ A character string giving the title of the plot. }
  \item{xlab}{ A character string naming the \eqn{X} axis. }
  \item{\dots}{ Additional arguments affecting the print or the plot produced. }
}
\note{When the object of class "strata" contains a certain stratum, \code{plot.strata} removes from the data the units in this stratum before generating the histogram.}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{strata.bh}}, \code{\link{strata.cumrootf}}, \code{\link{strata.geo}}, \code{\link{strata.LH}} }
\examples{
data(MRTS)
cumrootf<-strata.cumrootf(x=MRTS,nclass=500,CV=0.01,Ls=4,alloc=c(0.5,0,0.5))
print(cumrootf)
plot(cumrootf)
plot(cumrootf,drop=5)
plot(cumrootf,logscale=TRUE)
geo<-strata.geo(x=MRTS,CV=0.01,Ls=4,alloc=c(0.5,0,0.5))
print(geo)
plot(geo,logscale=TRUE)
# The geometric method does not perform well because of some small units
LH<-strata.LH(x=MRTS,CV=0.01,Ls=4,alloc=c(0.5,0,0.5),takeall=1)
print(LH)
plot(LH,logscale=TRUE)
}
\keyword{ methods }
