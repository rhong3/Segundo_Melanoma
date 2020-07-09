# validate OLA/COX proteins and phosphosites using TCGA RNAseq data
RNAseq <- read.csv("Data/TCGA_RNA/log2_RNAseq.csv")
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


# Aggregate Results
trans_data=read.csv("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/ola_data.csv")
trans_data = trans_data[,2:3]
trans=read.delim("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/3-yr-survival/dead_comparison_qvals.txt")
trans= trans[trans$significant == "True", ]
trans = trans[,1:2]
trt = merge(trans, trans_data, by='Hugo_Symbol')
trt$Feature = "3-yr-survival"
trt$Enriched_in = "dead"
trt_all = trt

trans_data=read.csv("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/ola_data.csv")
trans_data = trans_data[,2:3]
trans=read.delim("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/1-yr-survival/dead_comparison_qvals.txt")
trans= trans[trans$significant == "True", ]
trans = trans[,1:2]
trt = merge(trans, trans_data, by='Hugo_Symbol')
trt$Feature = "1-yr-survival"
trt$Enriched_in = "dead"
trt_all = rbind(trt_all, trt)

trans_data=read.csv("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/ola_data.csv")
trans_data = trans_data[,2:3]
trans=read.delim("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/dead_comparison_qvals.txt")
trans= trans[trans$significant == "True", ]
trans = trans[,1:2]
trt = merge(trans, trans_data, by='Hugo_Symbol')
trt$Feature = "5-yr-survival"
trt$Enriched_in = "dead"
trt_all = rbind(trt_all, trt)

trans_data=read.csv("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/ola_data.csv")
trans_data = trans_data[,2:3]
trans=read.delim("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/5-yr-survival/alive_comparison_qvals.txt")
trans= trans[trans$significant == "True", ]
trans = trans[,1:2]
trt = merge(trans, trans_data, by='Hugo_Symbol')
trt$Feature = "5-yr-survival"
trt$Enriched_in = "alive"
trt_all = rbind(trt_all, trt)

write.csv(trt_all, '~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA_summary.csv', row.names = FALSE)


# t-test alive/dead
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
toplist_prot = c('AIMP1', 'CDK4', 'DDX11', 'PAEP', 'CTNND1', 'GPR126', 'PIK3CB', 'TEX30', 'IARS', 'MMP12', 
                 'OXSR1', 'PDZD11', 'NCS1', 'RIOK1', 'THADA', 'FADD', 'NTPCR', 'TTYH3', 'XYLB')
toplist_phos = c('ADAM10', 'FGA', 'MARCKS', 'AKAP12', 'BCKDK', 'SRRM1', 'HMOX1', 'NIBAN2', 'PRKAB1', 'SYMPK')
svv = c('5-yr-survival', '3-yr-survival', '1-yr-survival', '6-month-survival')

for (sett in svv){
  dat = read.csv(paste("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/", sett, "/ola_data.csv", sep = ''))
  alive =  read_csv(paste("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/", sett,"/G1.txt", sep = ''), col_names = FALSE)$X1
  dead = read_csv(paste("~/documents/Segundo_Melanoma/Results/TCGA_RNA/OLA/", sett,"/G2.txt", sep = ''), col_names = FALSE)$X1
  for (gp in c('prot', 'phos')){
    sub.dat.long = data.frame()
    if (gp == 'prot'){
      sub.dat = dat[dat$Hugo_Symbol %in% toplist_prot, ]
      row.names(sub.dat) = sub.dat$Hugo_Symbol
      sub.dat = data.frame(t(sub.dat))
      sub.dat = sub.dat[-c(1:3), ]
      sub.dat$survival[which(row.names(sub.dat) %in% alive)] = 'alive'
      sub.dat$survival[which(row.names(sub.dat) %in% dead)] = 'dead'
    }
    else{
      sub.dat = dat[dat$Hugo_Symbol %in% toplist_phos, ]
      row.names(sub.dat) = sub.dat$Hugo_Symbol
      sub.dat = data.frame(t(sub.dat))
      sub.dat = sub.dat[-c(1:3), ]
      sub.dat$survival[which(row.names(sub.dat) %in% alive)] = 'alive'
      sub.dat$survival[which(row.names(sub.dat) %in% dead)] = 'dead'
    }

    for (i in 2:ncol(sub.dat)-1){
      selected = sub.dat[, c(i, ncol(sub.dat))]
      selected$gene = colnames(sub.dat)[i]
      colnames(selected) = gsub(colnames(sub.dat)[i], 'expression', colnames(selected))
      sub.dat.long = rbind(sub.dat.long, selected)
    }  
    sub.dat.long[, 1] <- sapply(sub.dat.long[, 1], as.character)
    sub.dat.long[, 1] <- sapply(sub.dat.long[, 1], as.numeric)
    
    pp = ggboxplot(sub.dat.long, x = "gene", y = "expression",
                   color = "black", fill = "survival", palette = "grey", title = paste(gp, 'top genes in TCGA RNAseq', sett, sep = ' ')) +
      stat_compare_means(method.args = list(alternative = "two.sided"), aes(group = survival), label = "p.signif", label.y = 18) + 
      stat_compare_means(method.args = list(alternative = "two.sided"), aes(group = survival), label = "p.format", label.y = 19)+ 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    
    pdf(file=paste("~/documents/Segundo_Melanoma/Results/TCGA_RNA/", gp, sett, ".pdf", sep=''),
        width=20,height=5)
    grid.arrange(pp,nrow=1, ncol=1)
    dev.off()
  }
}




