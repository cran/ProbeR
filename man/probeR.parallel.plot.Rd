\name{probeR.parallel.plot}
\alias{probeR.parallel.plot}
\title{ Parallel coordinate plot of summary value and probe level data}
\description{
  Parallel coordinate plot of summary value and probe level data is very useful to explore
the probe level data with summary value gene by gene.
}
\usage{
probeR.parallel.plot(affy.ID,data.summary, data.probe)
}
\arguments{
  \item{affy.ID}{affymetrix id}
  \item{data.summary}{ the summary values from the function exprs }
  \item{data.probe}{the normalized probe level data. They should be normalized using the same
method for the summary values. Also they should be the same probe level data calculating
the summary values. If the user use RMA with PM values,  this should be only PM values.}
}
\references{ Using Reliability with Gene Expression Models }
\author{Eun-Kyung Lee, Dianne Cook, Heike Hofmann, Maneesha Aluru, and Steve Rodermel }

\seealso{  \code{\link{probeR}}, ~~~ }
\examples{
library(affy)
data(affybatch.example)
eset<-expresso(affybatch.example,bg.correct=FALSE,normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
data.summary<-exprs(eset)
probe.norm<-normalize.AffyBatch.quantiles(affybatch.example,type="pmonly")
data.probe<-log(probes(probe.norm),base=2)
probeR.parallel.plot("D78156_at",data.summary,data.probe)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% at least one, from doc/KEYWORDS

