## Summarize findings (ICA-GSEA & OLA)
#' 
#' Created on 1/15/2020
#' 
#' @author: RH


library(dplyr)
# Aggregate GSEA
Prot_GSEA = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/ICA/significant_IC_clinical.csv")
Trans_GSEA = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/significant_IC_clinical.csv")
phospho_GSEA = read.csv("~/documents/Segundo_Melanoma/Results/phospho/ICA/significant_IC_clinical.csv")
GSEA = data.frame(matrix(ncol=11, nrow=0))
colnames(GSEA) = c('pathway',	'pval',	'padj',	'ES',	'NES',	'nMoreExtreme',	'size',	'leadingEdge', 'group', 'IC', 'clinical')
for (a in 1:nrow(Prot_GSEA)){
  m = toString(droplevels(Prot_GSEA[a, "Feature"]))
  n = toString(droplevels(Prot_GSEA[a, "IC"]))
  ICS = read.csv(paste("~/documents/Segundo_Melanoma/Results/proteomics/GSEA/", n, ".csv", sep=''), row.names=1)
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
  ICS = read.csv(paste("~/documents/Segundo_Melanoma/Results/transcriptomics/GSEA/", n, ".csv", sep=''), row.names=1)
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
  ICS = read.csv(paste("~/documents/Segundo_Melanoma/Results/phospho/GSEA/", n, ".csv", sep=''), row.names=1)
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
write.csv(GSEA, file = "~/documents/Segundo_Melanoma/Results/ICA_GSEA_summary.csv", row.names=FALSE)

GSEA.p = GSEA[GSEA$group == "proteomics", ]
GSEA.t = GSEA[GSEA$group == "transcriptomics", ]
GSEA.h = GSEA[GSEA$group == "phospho", ]
colnames(GSEA.h)[2:10] = paste(colnames(GSEA.h)[2:10], '.phos', sep='')
GSEA.joined = merge(GSEA.p, GSEA.t, by=c("pathway","clinical"), suffixes = c(".prot",".trans"))
GSEA.joined = merge(GSEA.joined, GSEA.h, by=c("pathway","clinical"))
write.csv(GSEA.joined, file = "~/documents/Segundo_Melanoma/Results/ICA_GSEA_summary_joined.csv", row.names=FALSE)

# Aggregate OLA
prot_data=read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/ola_data.csv")
prot_data=prot_data[,2:4]
# trans_data=read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/ola_data.csv")
# trans_data = trans_data[,2:3]
phospho_data=read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/ola_data.csv")
phospho_data = phospho_data[,2:5]
prot=read.delim("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5yr-survival-compare_alive_comparison_qvals.txt")
prot= prot[prot$significant == "True", ]
prot = prot[,1:2]
# trans=read.delim("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5yr-survival-compare_alive_comparison_qvals.txt")
# trans= trans[trans$significant == "True", ]
# trans = trans[,1:2]
phospho=read.delim("~/documents/Segundo_Melanoma/Results/phospho/OLA/5yr-survival-compare_alive_comparison_qvals.txt")
phospho= phospho[phospho$significant == "True", ]
phospho = phospho[,1:2]
ppt = merge(prot, prot_data, by="Accession")
ppt$Group = "proteomics"
# trt = merge(trans, trans_data, by="Gene.name")
# trt$Group = "transcriptomics"
pht = merge(phospho, phospho_data, by="Modified_sequence")
pht$Group = "phospho"
ddt = rbind.fill(ppt, pht)
ddt$Feature = "5yr-survival"
ddt$Enriched_in = "alive"
ddt = ddt[,c(4,3,1,5,2,6,7)]
ddt = ddt[order(ddt$FDR), ]
write.csv(ddt, file = "~/documents/Segundo_Melanoma/Results/OLA_summary.csv", row.names=FALSE)



