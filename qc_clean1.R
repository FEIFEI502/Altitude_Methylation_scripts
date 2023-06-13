library(meffil)
library(dplyr)
library(data.table)
options(mc.cores = 16)
setwd("##")
samplesheet=read.csv("##/samplesheet_qc_step4.csv") %>%
mutate(ori_name=Sample_Name) %>% mutate(Sample_Name=Sample_ID)
samplesheet$Slide = as.factor(samplesheet$Slide)
samplesheet$sentrix_row = as.factor(samplesheet$sentrix_row)
samplesheet$sentrix_col = as.factor(samplesheet$sentrix_col)
genotypes <- meffil.extract.genotypes("##/genotypes_QC.raw")
genotypes <- genotypes[,match(samplesheet$Sample_Name, colnames(genotypes))]

qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=TRUE)

##1802 samples
save(qc.objects,file="##/qc.objects.Robj")

qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold = 0.1,
  detectionp.samples.threshold = 0.1,
  detectionp.cpgs.threshold = 0.1,
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 5,
	snp.concordance.threshold = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters
)
save(qc.summary, file="##/qcsummary_c1.Robj")
meffil.qc.report(qc.summary, output.file="##/qc-report_c1.html")
