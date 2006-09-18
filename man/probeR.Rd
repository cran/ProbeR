\name{probeR}
\alias{probeR}
\title{ Calculate the reliability of one gene }
\description{
  calculate the reliability of one gene using summary variance and probe data variance
}
\usage{
probeR(summary.value, probe.value)
}
\arguments{
  \item{summary.value}{ the summary values of a specific gene for all the chips}
  \item{probe.value}{the probe level data of a specific gene. It is a matrix that the number of row is the number of
probes and the number of column is the number of chips}
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

\seealso{  \code{\link{probeR.wholegene}}, ~~~ }
\examples{
library(affy)
data(affybatch.example)
eset<-expresso(affybatch.example,bg.correct=FALSE,normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
data.summary<-exprs(eset)
probe.norm<-normalize.AffyBatch.quantiles(affybatch.example,type="pmonly")
data.probe<-log(probes(probe.norm),base=2)
summary.value<-data.summary[1,]
probe.value<-data.probe[1:16,]
probeR(summary.value,probe.value)
data.rel<-probeR.wholegene(data.summary,data.probe)

}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% at least one, from doc/KEYWORDS

