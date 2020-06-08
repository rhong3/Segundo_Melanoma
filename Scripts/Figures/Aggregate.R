## Summarize findings (ICA-GSEA & OLA)
#' 
#' Created on 1/15/2020
#' 
#' @author: RH


library(dplyr)
# Aggregate GSEA
Prot_GSEA = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/ICA/0.00001/significant_IC_clinical.csv")
Trans_GSEA = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/0.00001/significant_IC_clinical.csv")
phospho_GSEA = read.csv("~/documents/Segundo_Melanoma/Results/phospho/ICA/0.00001/significant_IC_clinical.csv")
GSEA = data.frame(matrix(ncol=11, nrow=0))
colnames(GSEA) = c('pathway',	'pval',	'padj',	'ES',	'NES',	'nMoreExtreme',	'size',	'leadingEdge', 'group', 'IC', 'clinical')
for (a in 1:nrow(Prot_GSEA)){
  m = toString(droplevels(Prot_GSEA[a, "Feature"]))
  n = toString(droplevels(Prot_GSEA[a, "IC"]))
  ICS = read.csv(paste("~/documents/Segundo_Melanoma/Results/proteomics/GSEA/strict/", n, ".csv", sep=''), row.names=1)
  ICS = ICS[ICS$padj < 0.01,]
  if (nrow(ICS) != 0){
    ICS$group = 'proteomics'
    ICS$IC = n
    ICS$clinical = m
    GSEA = rbind(GSEA, ICS)
  }
}

for (a in 1:nrow(Trans_GSEA)){
  m = toString(droplevels(Trans_GSEA[a, "Feature"]))
  n = toString(droplevels(Trans_GSEA[a, "IC"]))
  ICS = read.csv(paste("~/documents/Segundo_Melanoma/Results/transcriptomics/GSEA/strict/", n, ".csv", sep=''), row.names=1)
  ICS = ICS[ICS$padj < 0.01,]
  if (nrow(ICS) != 0){
    ICS$group = 'transcriptomics'
    ICS$IC = n
    ICS$clinical = m
    GSEA = rbind(GSEA, ICS)
  }
}

for (a in 1:nrow(phospho_GSEA)){
  m = toString(droplevels(phospho_GSEA[a, "Feature"]))
  n = toString(droplevels(phospho_GSEA[a, "IC"]))
  ICS = read.csv(paste("~/documents/Segundo_Melanoma/Results/phospho/GSEA/strict/", n, ".csv", sep=''), row.names=1)
  ICS = ICS[ICS$padj < 0.01,]
  if (nrow(ICS) != 0){
    ICS$group = 'phospho'
    ICS$IC = n
    ICS$clinical = m
    GSEA = rbind(GSEA, ICS)
  }
}
GSEA$clinical = as.character(GSEA$clinical)
GSEA = GSEA[order(GSEA$padj), ]
write.csv(GSEA, file = "~/documents/Segundo_Melanoma/Results/strict_ICA_GSEA_summary.csv", row.names=FALSE)

GSEA.p = GSEA[GSEA$group == "proteomics", ]
GSEA.t = GSEA[GSEA$group == "transcriptomics", ]
GSEA.h = GSEA[GSEA$group == "phospho", ]
colnames(GSEA.h)[2:10] = paste(colnames(GSEA.h)[2:10], '.phos', sep='')
GSEA.joined = merge(GSEA.p, GSEA.t, by=c("pathway","clinical"), suffixes = c(".prot",".trans"))
GSEA.joined = merge(GSEA.joined, GSEA.h, by=c("pathway","clinical"))
write.csv(GSEA.joined, file = "~/documents/Segundo_Melanoma/Results/strict_ICA_GSEA_summary_joined.csv", row.names=FALSE)


library(plyr)
# Aggregate OLA
# 5yr survival
# prot_data=read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/ola_data.csv")
# prot_data=prot_data[,2:4]
# trans_data=read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/ola_data.csv")
# trans_data = trans_data[,2:3]
# phospho_data=read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/ola_data.csv")
# phospho_data = phospho_data[,2:5]
# prot=read.delim("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/5-yr-survival-compare_alive_comparison_qvals.txt")
# prot= prot[prot$significant == "True", ]
# prot = prot[,1:2]
# trans=read.delim("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5yr-survival-compare_alive_comparison_qvals.txt")
# trans= trans[trans$significant == "True", ]
# trans = trans[,1:2]
# phospho=read.delim("~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/5-yr-survival-compare_alive_comparison_qvals.txt")
# phospho= phospho[phospho$significant == "True", ]
# phospho = phospho[,1:2]
# ppt = merge(prot, prot_data, by="Accession")
# ppt$Group = "proteomics"
# trt = merge(trans, trans_data, by="Gene.name")
# trt$Group = "transcriptomics"
# pht = merge(phospho, phospho_data, by="Modified_sequence")
# pht$Group = "phospho"
# ddt = rbind.fill(ppt, pht)
# ddt$Feature = "5-yr-survival"
# ddt$Enriched_in = "alive"
# ddt = ddt[,c(4,3,1,5,2,6,7,8)]
# ddt = ddt[order(ddt$FDR), ]

# 6 month survival
prot_data=read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/ola_data.csv")
prot_data=prot_data[,2:4]
trans_data=read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/ola_data.csv")
trans_data = trans_data[,2:3]
phospho_data=read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/ola_data.csv")
phospho_data = phospho_data[,2:5]
prot=read.delim("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/dead_comparison_qvals.txt")
# prot= prot[prot$significant == "True", ]
prot = prot[,1:2]
trans=read.delim("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/dead_comparison_qvals.txt")
# trans= trans[trans$significant == "True", ]
trans = trans[,1:2]
phospho=read.delim("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/dead_comparison_qvals.txt")
# phospho= phospho[phospho$significant == "True", ]
phospho = phospho[,1:2]
ppt = merge(prot, prot_data, by="Accession")
ppt$Group = "proteomics"
trt = merge(trans, trans_data, by="Gene.name")
trt$Group = "transcriptomics"
pht = merge(phospho, phospho_data, by="Modified_sequence")
pht$Group = "phospho"

ddt.2 = rbind.fill(ppt, pht, trt)
ddt.2$Feature = "6-month-survival"
ddt.2$Enriched_in = "dead"
ddt.2 = ddt.2[,c(4,3,1,5,2,6,7,8)]
ddt.2 = ddt.2[order(ddt.2$FDR), ]

ddt = ddt.2

# 1 year survival
prot_data=read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/ola_data.csv")
prot_data=prot_data[,2:4]
trans_data=read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/ola_data.csv")
trans_data = trans_data[,2:3]
phospho_data=read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/ola_data.csv")
phospho_data = phospho_data[,2:5]
prot=read.delim("~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/dead_comparison_qvals.txt")
# prot= prot[prot$significant == "True", ]
prot = prot[,1:2]
trans=read.delim("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/dead_comparison_qvals.txt")
# trans= trans[trans$significant == "True", ]
trans = trans[,1:2]
phospho=read.delim("~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/dead_comparison_qvals.txt")
# phospho= phospho[phospho$significant == "True", ]
phospho = phospho[,1:2]
ppt = merge(prot, prot_data, by="Accession")
ppt$Group = "proteomics"
trt = merge(trans, trans_data, by="Gene.name")
trt$Group = "transcriptomics"
pht = merge(phospho, phospho_data, by="Modified_sequence")
pht$Group = "phospho"

ddt.2 = rbind.fill(ppt, pht, trt)
ddt.2$Feature = "1-yr-survival"
ddt.2$Enriched_in = "dead"
ddt.2 = ddt.2[,c(4,3,1,5,2,6,7,8)]
ddt.2 = ddt.2[order(ddt.2$FDR), ]

ddt = rbind.fill(ddt, ddt.2)

# NRAS
prot_data=read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/ola_data.csv")
prot_data=prot_data[,2:4]
prot=read.delim("~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/mut_comparison_qvals.txt")
# prot= prot[prot$significant == "True", ]
prot = prot[,1:2]
ppt = merge(prot, prot_data, by="Accession")
ppt$Group = "proteomics"
ppt$Feature = "NRAS"
ppt$Enriched_in = "mut"

ddt = rbind.fill(ddt, ppt)

# stage
prot_data=read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/ola_data.csv")
prot_data=prot_data[,2:4]
prot=read.delim("~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/4_comparison_qvals.txt")
# prot= prot[prot$significant == "True", ]
prot = prot[,1:2]
ppt = merge(prot, prot_data, by="Accession")
ppt$Group = "proteomics"
ppt$Feature = "stage"
ppt$Enriched_in = "4"

ddt = rbind.fill(ddt, ppt)

write.csv(ddt, file = "~/documents/Segundo_Melanoma/Results/full_OLA_summary.csv", row.names=FALSE)

old=read.csv("~/documents/Segundo_Melanoma/Legacy/old_summary/OLA_summary_old.csv")

mgg = merge(x=ddt, y=old, by = c('Gene.name', 'Description', 'Group', 'Feature', 'Accession', 'Modified_sequence'), all = TRUE, suffixes=c('_new', '_old'))
write.csv(mgg, file="~/documents/Segundo_Melanoma/Legacy/old_summary/OLA_new_old.csv", row.names=FALSE)




