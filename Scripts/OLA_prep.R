## Outlier Analysis prep
#' 
#' Created on 1/15/2020
#' 
#' @author: RH

# proteomics
prot = read.csv("~/documents/Lund_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Lund_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[prot$os.events == 1 | prot$os.days >= 1825,]
prot.ola = prot.ola[!is.na(prot.ola$os.days),]
prot.ola = prot.ola[!is.na(prot.ola$os.events),]
write.csv(prot.ola, "~/documents/Lund_Melanoma/Results/proteomics/OLA/5yr-survival-clinical.csv")
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Lund_Melanoma/Results/proteomics/OLA/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Lund_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Lund_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[trans$os.events == 1 | trans$os.days >= 1825,]
trans.ola = trans.ola[!is.na(trans.ola$os.days),]
trans.ola = trans.ola[!is.na(trans.ola$os.events),]
write.csv(trans.ola, "~/documents/Lund_Melanoma/Results/transcriptomics/OLA/5yr-survival-clinical.csv")
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Lund_Melanoma/Results/transcriptomics/OLA/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Lund_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Lund_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[phospho$os.events == 1 | phospho$os.days >= 1825,]
phospho.ola = phospho.ola[!is.na(phospho.ola$os.days),]
phospho.ola = phospho.ola[!is.na(phospho.ola$os.events),]
write.csv(phospho.ola, "~/documents/Lund_Melanoma/Results/phospho/OLA/5yr-survival-clinical.csv")
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Lund_Melanoma/Results/phospho/OLA/ola_data.csv")

