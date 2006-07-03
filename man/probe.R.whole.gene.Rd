\name{probe.R.whole.gene}
\alias{probe.R.whole.gene}
\title{ Calculate the reliability of the whole genes }
\description{
  calculate the reliability of the whole genes using summary variance and probe data variance
}
\usage{
probe.R.whole.gene(data.summary, data.probe)
}
\arguments{
  \item{data.summary}{ the summary values from the function exprs }
  \item{data.probe}{the normalized probe level data. They should be normalized using the same
method for the summary values. Also they should be the same probe level data calculating
the summary values. If the user use RMA with PM values,  this should be only PM values.}
}
\value{
  \item{summary.var }{ the variance of the summary values}
  \item{probe.var }{ the variance of the probe level data}
  \item{probe.n }{ the number of probes}
  \item{reliability }{ the reliability }
  ...
}
\references{ Using Reliability with Gene Expression Models }
\author{Eun-Kyung Lee, Dianne Cook, Heike Hofmann, Maneesha Aluru, and Steve Rodermel }

\seealso{  \code{\link{probe.R}}, ~~~ }
\examples{
library(affy)
data(affybatch.example)
eset<-expresso(affybatch.example,bg.correct=FALSE,normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
data.summary<-exprs(eset)
probe.norm<-normalize.AffyBatch.quantiles(affybatch.example,type="pmonly")
data.probe<-log(probes(probe.norm),base=2)
summary.value<-data.summary[1,]
probe.value<-data.probe[1:16,]
probe.R(summary.value,probe.value)
data.rel<-probe.R.whole.gene(data.summary,data.probe)

}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% at least one, from doc/KEYWORDS

