library(meffil)
library(dplyr)
options(mc.cores = 16)

load("qc.objects.c1.Robj")
load("qcsummary.c2.Robj")
samplesheet=read.csv("samplesheet_qc_step5.csv")
samplesheet$Slide = as.factor(samplesheet$Slide)
samplesheet$sentrix_row = as.factor(samplesheet$sentrix_row)
samplesheet$sentrix_col = as.factor(samplesheet$sentrix_col)

pc = 20  
norm.objects = meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste("##/norm.obj.pc",pc,".Robj",sep=""))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta,file=paste("##/norm.beta.pc",pc,".Robj",sep=""))

cc = t(sapply(qc.objects, function(obj)
  obj$cell.counts$counts))
cc = data.frame(IID=row.names(cc),cc)
write.table(cc,"##/cellcounts_combined.txt",sep="\t", row.names=F,col.names=T,quote=F)                                                         



batch_var<-c("Slide",
             "sentrix_row",
             "sentrix_col",
             "Sex","batch")
norm.parameters <- meffil.normalization.parameters( norm.objects,
                                                    variables=batch_var,
                                                    control.pcs=1:10,
                                                    batch.pcs=1:10,
                                                    batch.threshold=0.01
)

pcs = meffil.methylation.pcs(norm.beta, probe.range=20000)
save(pcs,file="##/pcs.norm.beta_combined.Robj")
norm.summary= meffil.normalization.summary(norm.objects, pcs=pcs,parameters=norm.parameters) 
meffil.normalization.report(norm.summary, output.file="normalization-report_combined.html")
