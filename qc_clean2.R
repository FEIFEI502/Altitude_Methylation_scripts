library(meffil)
library(dplyr)
library(data.table)

load("qc.objects.Robj")
load("qcsummary_c1.Robj")
outlier = qc.summary$bad.samples
write.table(outlier, file="outlier_c1.txt", quote=F, sep="\t")

length(qc.objects)
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects)

save(qc.objects,file="qc.objects.c1.Robj")

qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold = 0.1,
  detectionp.samples.threshold = 0.1,
  detectionp.cpgs.threshold = 0.1,
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 5,
	snp.concordance.threshold = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

samplesheet=read.csv("samplesheet_qc_step5.csv") 
samplesheet$Slide = as.factor(samplesheet$Slide)
samplesheet$sentrix_row = as.factor(samplesheet$sentrix_row)
samplesheet$sentrix_col = as.factor(samplesheet$sentrix_col)
genotypes <- meffil.extract.genotypes("##/genotypes_QC.raw")
genotypes <- genotypes[,match(samplesheet$Sample_Name, colnames(genotypes))]


qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters,
genotypes=genotypes)

save(qc.summary, file="qcsummary.c2.Robj")
meffil.qc.report(qc.summary,output.file="qc-report.c2.html")
