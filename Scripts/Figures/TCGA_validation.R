# validate OLA/COX proteins and phosphosites using TCGA RNAseq data
RNAseq <- read_csv("Data/TCGA_RNA/log2_RNAseq.csv")
clinical <- read.delim("~/Documents/Segundo_melanoma/Data/TCGA_RNA/data_clinical_sample.txt", comment.char="#")
clinical$SAMPLE_ID=gsub('-','.', clinical$SAMPLE_ID)
patient = read.delim("~/Documents/Segundo_melanoma/Data/TCGA_RNA/data_clinical_patient.txt", comment.char="#")
patient.clinical = merge(patient, clinical, by="PATIENT_ID")
row.names(patient.clinical) = patient.clinical$SAMPLE_ID

patient.clinical = patient.clinical[, c('DSS_MONTHS', 'DSS_STATUS')]
patient.clinical = na.omit(patient.clinical)
patient.clinical$survival_6mo = NA
patient.clinical$survival_6mo[which(patient.clinical$DSS_MONTHS >= 6)] = 1
patient.clinical$survival_6mo[which(patient.clinical$DSS_MONTHS < 6 & patient.clinical$DSS_STATUS == 'DEAD WITH TUMOR')] = 0
patient.clinical$survival_1yr=NA
patient.clinical$survival_1yr[which(patient.clinical$DSS_MONTHS >= 12)] = 1
patient.clinical$survival_1yr[which(patient.clinical$DSS_MONTHS < 12 & patient.clinical$DSS_STATUS == 'DEAD WITH TUMOR')] = 0
patient.clinical$survival_3yr=NA
patient.clinical$survival_3yr[which(patient.clinical$DSS_MONTHS >= 36)] = 1
patient.clinical$survival_3yr[which(patient.clinical$DSS_MONTHS < 36 & patient.clinical$DSS_STATUS == 'DEAD WITH TUMOR')] = 0
patient.clinical$survival_5yr=NA
patient.clinical$survival_5yr[which(patient.clinical$DSS_MONTHS >= 60)] = 1
patient.clinical$survival_5yr[which(patient.clinical$DSS_MONTHS < 60 & patient.clinical$DSS_STATUS == 'DEAD WITH TUMOR')] = 0

RNAseq_survival = RNAseq[, colnames(RNAseq) %in% c(rownames(patient.clinical), 'Hugo_Symbol', 'Entrez_Gene_Id')]
RNAseq_survival = na.omit(RNAseq_survival)

trans.ola = patient.clinical[!is.na(patient.clinical$survival_5yr),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/5yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_5yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_5yr==0),]), fileConn)
close(fileConn)
trans.oladata = RNAseq_survival[, intersect(rownames(trans.ola), colnames(RNAseq_survival))]
trans.oladata = cbind(RNAseq_survival[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/ola_data.csv")


trans.ola = patient.clinical[!is.na(patient.clinical$survival_3yr),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/3yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_3yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_3yr==0),]), fileConn)
close(fileConn)
trans.oladata = RNAseq_survival[, intersect(rownames(trans.ola), colnames(RNAseq_survival))]
trans.oladata = cbind(RNAseq_survival[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/ola_data.csv")


trans.ola = patient.clinical[!is.na(patient.clinical$survival_1yr),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/1yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_1yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_1yr==0),]), fileConn)
close(fileConn)
trans.oladata = RNAseq_survival[, intersect(rownames(trans.ola), colnames(RNAseq_survival))]
trans.oladata = cbind(RNAseq_survival[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/ola_data.csv")


trans.ola = patient.clinical[!is.na(patient.clinical$survival_6mo),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/6-month-survival/6mo-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/6-month-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/6-month-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_6mo==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/6-month-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_6mo==0),]), fileConn)
close(fileConn)
trans.oladata = RNAseq_survival[, intersect(rownames(trans.ola), colnames(RNAseq_survival))]
trans.oladata = cbind(RNAseq_survival[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/6-month-survival/ola_data.csv")



