##MWAS in batch

library(data.table)
library(dplyr)
#library(PCAForQTL)
args = commandArgs(trailingOnly=TRUE)

##input raw meth data
setwd("/storage/yangjianLab/chengfeifei/130_tibetan_methylation/130_tibetan_meth_v3/meth")
load("norm.pc20.autosome.beta.Robj") # meth 445001 * 1622

samplesheet=fread("samplesheet_qc_step6.csv")  %>%
          mutate(STM=case_when(Group=="TPHan"~1,Group=="WZHan"~0)) %>%
            mutate(LTM=case_when(Group=="TPTib"~1,Group=="TPHan"~0))

samplesheet$Slide = as.factor(samplesheet$Slide)
samplesheet$sentrix_row = as.factor(samplesheet$sentrix_row)
samplesheet$sentrix_col = as.factor(samplesheet$sentrix_col)

####select corresponding batch and samples
subgroup= samplesheet %>% filter(batch==as.numeric(args[1])) %>% filter(! is.na(args[2]))
meth_subgroup= subset(meth,select=subgroup$Sample_Name)

###check the sampleID matched
print("check whether the sampleID is matched")
table(colnames(meth_subgroup)==subgroup$Sample_Name)

##cell component
cc=fread("cellcounts_combined.txt") %>% filter(IID %in% subgroup$Sample_Name) 
cc_subgroup = cc %>% select(-IID)
cc_subgroup= as.matrix(cc_subgroup)

print("check whether the sampleID is matched")
table(cc$IID==subgroup$Sample_Name)

##methylation PCs
#mpc=prcomp(t(meth_subgroup),scale=T)
#PCs<-mpc$x
#resultRunElbow<-PCAForQTL::runElbow(prcompResult=mpc)
#print(paste("resultRunElbow is", resultRunElbow,sep=" ")

#PCsTop<-as.matrix(PCs[,1:resultRunElbow])

####perform EWAS
n = nrow(meth_subgroup)
sto = matrix(NA, nrow = n, ncol = 4)

for(i in 1:n) {
if(i%%100 == 0) print(i)
fm <- as.formula(paste("meth_subgroup[i,] ~", args[2],"+ age + Sex + Slide + sentrix_row + sentrix_col + cc_subgroup"))
temp = lm(fm,data=subgroup)
sto[i,] = summary(temp)$coef[2,]
}

result= as.data.frame(sto)
colnames(result)=c("beta","SE","t","p")
result$Probe= rownames(meth)

write.table(result,paste("/storage/yangjianLab/chengfeifei/130_tibetan_methylation/130_tibetan_meth_v3/EWAS/batch_",args[1],"/",args[2],"_EWAS_basicNoMethPC.txt",sep=""),,quote=F,row.names=F,col.names=T)
