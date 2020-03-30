library(dplyr)
library(plyr)

OLA = read.csv("~/Documents/Segundo_Melanoma/Results/OLA_summary.csv")
prot.clinical = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names = 1)
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
phospho.clinical = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names = 1)
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")

OLA.1 = OLA[OLA['Group'] == "phospho",]
OLA.1=OLA.1[OLA.1['Feature'] == '5-yr-survival', ]
prot.1 = subset(prot, Accession %in% c('P22059', 'Q9Y2I1', 'Q8IWB9', 'Q9UQ35'))
phospho.1 = subset(phospho, Modified_sequence %in% c('_MLAES[Phospho (STY)]DESGDEESVSQTDK_', '_TPS[Phospho (STY)]PEPVDK_', '_TAPSS[Phospho (STY)]PLTSPSDTR_', '_SLS[Phospho (STY)]YSPVER_'))
phospho.clinical$name=rownames(phospho.clinical)
prot.clinical$name=rownames(prot.clinical)

clinical = na.omit(merge(phospho.clinical, prot.clinical)[, c(43,44)])
nn = intersect(clinical$name, colnames(prot.1))
prot.2 = prot.1[, c(nn, 'Accession', 'Gene.name')]
mm = intersect(clinical$name, colnames(phospho.1))
phospho.2 = phospho.1[, c(mm, 'Accession', 'Gene.name', 'Modified_sequence')]
rownames(phospho.2) = phospho.2$Accession
rownames(phospho.2) = paste(rownames(phospho.2), 'phos', sep='_')
rownames(prot.2) = prot.2$Accession
rownames(prot.2) = paste(rownames(prot.2), 'prot', sep='_')
phospho.3 = t(phospho.2)
prot.3 = t(prot.2)
phospho.3 = phospho.3[-c(67,68,69),]
prot.3 = prot.3[-c(67,68),]
joined = cbind(prot.3, phospho.3)
rownames(clinical) = clinical$name
joined = merge(joined, clinical, by=0)
rownames(joined) = joined$Row.names
joined = joined[,-c(1,11)] 

fff = c('P22059_phos', 'Q9Y2I1_phos', 'Q9UQ35_phos', 'Q8IWB9_phos')
for (f in fff){
  print(f)
  x1 = as.numeric(as.character(unlist(joined[joined$survival_5yr == 0, ][f])))
  x2 = as.numeric(as.character(unlist(joined[joined$survival_5yr == 1, ][f])))
  print(t.test(x1, x2, alternative = "two.sided", var.equal = TRUE))
  x1a = x1-median(x1)
  x2a = x2-median(x2)
}

rownames(phospho.1) = phospho.1$Gene.name
rownames(prot.1) = prot.1$Gene.name
phospho.clinicalx = na.omit(phospho.clinical[,c(43,44)])
prot.clinicalx = na.omit(prot.clinical[,c(43,44)])
phosm = merge(data.frame(t(phospho.1)), data.frame(phospho.clinicalx), by=0, all=FALSE)[,1:6]
protm = merge(data.frame(t(prot.1)), data.frame(prot.clinicalx), by=0, all=FALSE)[,1:6]
protm = protm[, order(colnames(protm))]
phosm = phosm[, order(colnames(phosm))]

out = data.frame(matrix(nrow=4, ncol=4))
colnames(out) = c('prot_mean_Alive-Dead', 'prot_median_Alive-Dead', 'phos_mean_Alive-Dead', 'phos_median_Alive-Dead')
rownames(out) = c('NISCH', 'OSBP', 'SRRM2', 'TEX2')
for(m in c('NISCH', 'OSBP', 'SRRM2', 'TEX2')){
  median_alivep=median(as.numeric(as.character(protm[protm$survival_5yr==1,m])))
  median_deadp=median(as.numeric(as.character(protm[protm$survival_5yr==0,m])))
  medianfcp = median_alivep-median_deadp
  median_aliveh=median(as.numeric(as.character(phosm[phosm$survival_5yr==1,m])))
  median_deadh=median(as.numeric(as.character(phosm[phosm$survival_5yr==0,m])))
  medianfch = median_aliveh-median_deadh
  
  mean_alivep=mean(as.numeric(as.character(protm[protm$survival_5yr==1,m])))
  mean_deadp=mean(as.numeric(as.character(protm[protm$survival_5yr==0,m])))
  meanfcp = mean_alivep-mean_deadp
  mean_aliveh=mean(as.numeric(as.character(phosm[phosm$survival_5yr==1,m])))
  mean_deadh=mean(as.numeric(as.character(phosm[phosm$survival_5yr==0,m])))
  meanfch = mean_aliveh-mean_deadh
  
  out[m, 'prot_mean_Alive-Dead'] = meanfcp
  out[m, 'phos_mean_Alive-Dead'] = meanfch
  out[m, 'prot_median_Alive-Dead'] = medianfcp
  out[m, 'phos_median_Alive-Dead'] = medianfch
  
}

