\name{probeR.openggobi.with}
\alias{probeR.openggobi.with}
\title{ Explore summary and probe level data with the reliability - dynamic version }
\description{
  With GGobi, the user can explore summary and probe level data with the reliability dynamically. The user
also can add their own results from other statistical analysis.
}
\usage{
probeR.openggobi.with(data.summary, data.probe,selected.id,data.user.define)
}
\arguments{
  \item{data.summary}{ the summary values from the function exprs }
  \item{data.probe}{the normalized probe level data. They should be normalized using the same
method for the summary values. Also they should be the same probe level data calculating
the summary values. If the user use RMA with PM values,  this should be only PM values.}
  \item{selected.id}{the user selected id for statistical analysis result using other method}
  \item{data.user.define}{the user defined data for statistical analysis result using other method}
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
probeR.openggobi.with(data.summary,data.probe)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% at least one, from doc/KEYWORDS

