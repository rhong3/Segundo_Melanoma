## ML prep
#' 
#' Created on 1/17/2020
#' 
#' @author: RH

sig = read.csv("~/documents/Segundo_Melanoma/Results/OLA_summary.csv")
sig = sig[sig$Feature=="6-month-survival", ]
clinical.p = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/6-month-survival-clinical.csv")
clinical.h = read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/6-month-survival-clinical.csv")
clinical.t = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/6-month-survival-clinical.csv")

clinical = merge(clinical.p, clinical.h)
clinical = merge(clinical, clinical.t)

pdata = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/ola_data.csv")
tdata = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/ola_data.csv")
hdata = read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/ola_data.csv")
data.p = pdata[pdata$Accession %in% sig$Accession, ]
data.t = tdata[tdata$Gene.name %in% sig$Gene.name, ]
data.h = hdata[hdata$Modified_sequence %in% sig$Modified_sequence, ]
data.p$Gene.name = paste(data.p$Gene.name, ".pro", sep='')
rownames(data.p) = data.p$Gene.name
data.p = t(data.p)
data.p = data.p[rownames(data.p) %in% clinical$X, ]
data.t$Gene.name = paste(data.t$Gene.name, ".tra", sep='')
rownames(data.t) = data.t$Gene.name
data.t = t(data.t)
data.t = data.t[rownames(data.t) %in% clinical$X, ]
data.h$Gene.name = paste(data.h$Gene.name, ".pho", sep='')
rownames(data.h) = data.h$Gene.name
data.h = t(data.h)
data.h = data.h[rownames(data.h) %in% clinical$X, ]
rownames(clinical) = clinical$X
dataa = cbind(data.p, data.t, data.h, clinical)
dataa = dataa[,-c(122:165)]

write.csv(dataa, "~/documents/Segundo_Melanoma/Results/ML/6-month-survival.csv", row.names=TRUE)


# All data
sig = read.csv("~/documents/Segundo_Melanoma/Results/OLA_summary.csv")
sig = sig[sig$Feature=="1-year-survival", ]
clinical.p = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/1-yr-survival-clinical.csv")
clinical.h = read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/1-yr-survival-clinical.csv")
clinical.t = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/1-yr-survival-clinical.csv")
clinical = merge(clinical.p, clinical.h)
clinical = merge(clinical, clinical.t)

pdata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/ICA_proteomics.csv")
rownames(pdata) = pdata$X
tdata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/ICA_transcriptomics.csv")
rownames(tdata) = tdata$X
hdata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/ICA_phospho.csv")
rownames(hdata) = hdata$X
data.p = t(pdata)
data.p = data.p[rownames(data.p) %in% clinical$X, ]
data.t = t(tdata)
data.t = data.t[rownames(data.t) %in% clinical$X, ]
data.h = t(hdata)
data.h = data.h[rownames(data.h) %in% clinical$X, ]
rownames(clinical) = clinical$X
datap = cbind(data.p, clinical[45])
datat = cbind(data.t, clinical[45])
datah = cbind(data.h, clinical[45])

write.csv(datap, "~/documents/Segundo_Melanoma/Results/ML/1-yr-survival_proteomics.csv", row.names=TRUE)
write.csv(datat, "~/documents/Segundo_Melanoma/Results/ML/1-yr-survival_transcriptomics.csv", row.names=TRUE)
write.csv(datah, "~/documents/Segundo_Melanoma/Results/ML/1-yr-survival_phospho.csv", row.names=TRUE)

