probeR<-function(summary.value,probe.value)
{
summary.var<-var(summary.value)
probe.var<-mean(diag(var(probe.value)))
probe.n<-nrow(probe.value)
reliability<-summary.var/(summary.var+probe.var/probe.n)
return(list(summary.var=summary.var,probe.var=probe.var,
            probe.n=probe.n,reliability=reliability))
}


probeR.wholegene<-function(data.summary,data.probe)
{
#print(10)
   gene.names.probe<-matrix(rownames(data.probe))
   A<-matrix(apply(gene.names.probe,1,function(x){paste(strsplit(x,"_at")[[1]][1],"_at",sep="")}))
   ST.id<-which(is.na(apply(gene.names.probe,1,function(x){strsplit(x,"_at")[[1]][2]})))
   A[ST.id,]<-apply(as.matrix(gene.names.probe[ST.id,]),1,function(x){paste(strsplit(x,"_st")[[1]][1],"_st",sep="")})
   probe.length<-table(A)
    gene.names<-rownames(data.summary)


    rel.data<-NULL

    for(i in 1:length(probe.length))
    {  id.i<-i
       ifelse(id.i==1,id<-1,id<-sum(probe.length[1:(id.i-1)])+1)
       probe.value<-data.probe[id:(id+probe.length[id.i]-1),]
       summary.value<-data.summary[i,]
       rel.temp<-probeR(summary.value,probe.value)
       rel.data<-rbind(rel.data,c(rel.temp$summary.var,rel.temp$probe.var,
                        rel.temp$probe.n,rel.temp$reliability))
    }
    rownames(rel.data)<-gene.names
    colnames(rel.data)<-c("summary.var","probe.var","probe.n","reliability")
    return(rel.data)

}

probeR.openggobi.with<-function(data.summary,data.probe,selected.id=NA,data.user.define=NA)
{

   gene.names.PM<-matrix(rownames(data.probe))
   A<-matrix(apply(gene.names.PM,1,function(x){paste(strsplit(x,"_at")[[1]][1],"_at",sep="")}))
   ST.id<-which(is.na(apply(gene.names.PM,1,function(x){strsplit(x,"_at")[[1]][2]})))
   A[ST.id,]<-apply(as.matrix(gene.names.PM[ST.id,]),1,function(x){paste(strsplit(x,"_st")[[1]][1],"_st",sep="")})
   PM.length<-table(A)
   if(!is.na(selected.id))
   {  data.summary<-data.summary[selected.id,]
      new.data.probe<-c()
      probe.id<-c()
      for(i in 1:length(selected.id))
      {  new.data.probe<-rbind(new.data.probe,data.probe[A==rownames(data.summary)[i],])
         probe.id<-c(probe.id,rep(rownames(data.summary)[i],PM.length[selected.id[i]]))
      }
      data.probe<-new.data.probe
      data.user.define<-data.user.define[selected.id,]
   } else
   { probe.id<-c(A)
   }
   data.rel<-probeR.wholegene(data.summary,data.probe)
   if(is.na(data.user.define))
   { summary.ggobi<-data.frame(data.summary,reliability=data.rel,id=factor(rownames(data.summary)))
   } else
   { if(length(data.user.define)==nrow(data.summary) || nrow(data.user.define)==nrow(data.summary))
     { summary.ggobi<-data.frame(data.summary,reliability=data.rel,user=data.user.define,id=factor(rownames(data.summary)))
     } else
     { print(" Error : It should be length(data.user.define)==nrow(data.summary) or nrow(data.user.define)==nrow(data.summary)")
     }
   }
   probe.ggobi<-data.frame(data.probe,id=probe.id)

   require(rggobi)
 
   ggobi(summary.ggobi)

   ggobi_set_data_frame(probe.ggobi)

}

probeR.parallel.plot<-function(affy.ID,data.summary, data.probe)
{
   gene.names.PM<-matrix(rownames(data.probe))
   id<-which(rownames(data.summary)==affy.ID)
   A<-matrix(apply(gene.names.PM,1,function(x){paste(strsplit(x,"_at")[[1]][1],"_at",sep="")}))
   ST.id<-which(is.na(apply(gene.names.PM,1,function(x){strsplit(x,"_at")[[1]][2]})))
   A[ST.id,]<-apply(as.matrix(gene.names.PM[ST.id,]),1,function(x){paste(strsplit(x,"_st")[[1]][1],"_st",sep="")})
   probe.id<-which(A==affy.ID)
   xlab<-colnames(data.summary)
   probe.value<-data.probe[probe.id,]
   summary.value<-data.summary[id,]
   p<-ncol(data.summary)

   plot(1:p,summary.value,type='l',lty=1,lwd=2,ylim=range(data.probe),xlab="experiments",ylab="expression level",axes=FALSE)
   axis(2)
   axis(1,at=c(1:p),lab=as.character(xlab),cex=0.7)
   for(i in 1:length(probe.id))
     lines(c(1:p),probe.value[i,],lty=2)
   title(main=as.character(affy.ID))

} 
 .onLoad<-function(lib,pkg){
   require("rggobi",  quietly=TRUE) || stop("rggobi package not found")
}

