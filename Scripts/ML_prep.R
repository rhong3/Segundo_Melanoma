## ML prep
#' 
#' Created on 1/17/2020
#' 
#' @author: RH

sig = read.csv("~/documents/Segundo_Melanoma/Results/OLA_summary.csv")
clinical.p = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5yr-survival-clinical.csv")
clinical.p$survival = as.numeric(clinical.p$os.days>1825)
clinical.p = clinical.p[,c(1,44)]
clinical.h = read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/5yr-survival-clinical.csv")
clinical.h$survival = as.numeric(clinical.h$os.days>1825)
clinical.h = clinical.h[,c(1,44)]
clinical.t = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5yr-survival-clinical.csv")
clinical.t$survival = as.numeric(clinical.t$os.days>1825)
clinical.t = clinical.t[,c(1,44)]

clinical = merge(clinical.p, clinical.h)
clinical = merge(clinical, clinical.t)

pdata = read.csv("~/documents/Segundo_Melanoma/Results/proteomics/OLA/ola_data.csv")
tdata = read.csv("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/ola_data.csv")
hdata = read.csv("~/documents/Segundo_Melanoma/Results/phospho/OLA/ola_data.csv")
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
dataa = dataa[,-20]

write.csv(dataa, "~/documents/Segundo_Melanoma/Results/ML/5yr-survival.csv", row.names=TRUE)




