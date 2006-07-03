
probe.R<-function(summary.value, probe.value)
{
  summary.var<-var(summary.value)
  probe.var<-mean(diag(var(probe.value)))
  probe.n<-nrow(probe.value)
  reliability<-summary.var/(summary.var+probe.var/probe.n)
  return(list(summary.var=summary.var,probe.var=probe.var,
             probe.n=probe.n,reliability=reliability))
}


probe.R.whole.gene<-function(data.summary,data.probe)
{
    gene.names.probe<-matrix(rownames(data.probe))
    gene.names.probe<-apply(as.matrix(gene.names.probe),1,
                          function(x){strsplit(x,"_at")[[1]][1]})
    probe.length<-table(gene.names.probe)
    gene.names<-rownames(data.summary)
    rel.data<-NULL

    for(i in 1:length(probe.length))
    {  id.i<-i
       ifelse(id.i==1,id<-1,id<-sum(probe.length[1:(id.i-1)])+1)
       probe.value<-data.probe[id:(id+probe.length[id.i]-1),]
       summary.value<-data.summary[i,]
       rel.temp<-probe.R(summary.value,probe.value)
       rel.data<-rbind(rel.data,c(rel.temp$summary.var,rel.temp$probe.var,
                        rel.temp$probe.n,rel.temp$reliability))
    }
    rownames(rel.data)<-gene.names
    colnames(rel.data)<-c("summary.var","probe.var","probe.n","reliability")
    return(rel.data)
}
