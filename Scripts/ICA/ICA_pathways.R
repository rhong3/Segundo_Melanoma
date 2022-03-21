### ICA pathways ###
library(ggplot2)
library(tsne)
library(cluster)
library(reshape)

pathways <- read_excel("Data/20220321_ICA-GSEA_selected_pathways_v1.xlsx")
strict_ICA_GSEA_summary <- read.csv("~/Documents/Segundo_Melanoma/Results/strict_ICA_GSEA_summary.csv")
strict_ICA_GSEA_summary$logp = -log10(strict_ICA_GSEA_summary$padj)

prot = strict_ICA_GSEA_summary[strict_ICA_GSEA_summary$group=='proteomics' & strict_ICA_GSEA_summary$pathway %in% pathways$Proteomics, c('pathway', 'IC', 'logp')]
icPave = read.delim("~/Documents/Segundo_Melanoma/Results/Figures/Figure2-S2/ICA_proteomics_IC_Clinical_Correlation_P_Value_all.tsv", row.names = 1)
icPave = as.matrix(icPave)

icPave_plot=melt(icPave)
icPave_plot=subset(icPave_plot,
                   X1%in%names(which(apply(icPave,1,sum)>30))&X2%in%colnames(icPave)[which(apply(icPave,2,sum)>30)])
icPave_plot$value[icPave_plot$value > 100] = 100

colnames(icPave_plot)[1:2]=c('IC','clinical_feature')
icPave_plot[,1]=factor(icPave_plot[,1])
icPave_plot$IC=paste('IC',str_pad(icPave_plot$IC,3,pad='0'),sep='_')
IC = unique(icPave_plot$IC)[!unique(icPave_plot$IC) %in% prot$IC]
dum = as.data.frame(IC)
dum$pathway = NA
dum$logp = NA
dum = dum[, c("pathway", "IC", "logp")]
prot$IC = as.character(prot$IC)
prot = rbind(prot, dum)
prot = prot[order(prot$IC), ]
rownames(prot) = 1:nrow(prot)

exp_prefix='ICA_proteomics'
out_dir='~/documents/Segundo_Melanoma/Results/Figures'
png(filename=paste(out_dir,paste(exp_prefix,'IC_pathway.png',sep='_'),sep='/'),
    width = 1500, height = 500, units = "px", 
    pointsize = 20,
    bg = 'transparent')
ggplot(prot,aes(x=IC,y=pathway))+
  theme_classic(base_size=20)+
  theme(axis.line=element_blank())+
  theme(axis.text.x = element_text(angle=90,vjust=0.5, size=20))+
  theme(axis.text.y = element_text(vjust=0.5, size=20))+
  geom_tile(aes(fill=logp))+
  scale_fill_gradient2(low='white',high='purple')+
  geom_text(aes(label=round(logp,digits=2)),size=6)
graphics.off()



phospho = strict_ICA_GSEA_summary[strict_ICA_GSEA_summary$group=='phospho' & strict_ICA_GSEA_summary$pathway %in% pathways$Phospho, c('pathway', 'IC', 'logp')]
icPave = read.delim("~/Documents/Segundo_Melanoma/Results/Figures/Figure2-S2/ICA_phospho_IC_Clinical_Correlation_P_Value_all.tsv", row.names = 1)
icPave = as.matrix(icPave)

icPave_plot=melt(icPave)
icPave_plot=subset(icPave_plot,
                   X1%in%names(which(apply(icPave,1,sum)>30))&X2%in%colnames(icPave)[which(apply(icPave,2,sum)>30)])
icPave_plot$value[icPave_plot$value > 100] = 100

colnames(icPave_plot)[1:2]=c('IC','clinical_feature')
icPave_plot[,1]=factor(icPave_plot[,1])
icPave_plot$IC=paste('IC',str_pad(icPave_plot$IC,3,pad='0'),sep='_')
IC = unique(icPave_plot$IC)[!unique(icPave_plot$IC) %in% phospho$IC]
dum = as.data.frame(IC)
dum$pathway = NA
dum$logp = NA
dum = dum[, c("pathway", "IC", "logp")]
phospho$IC = as.character(phospho$IC)
phospho = rbind(phospho, dum)
phospho = phospho[order(phospho$IC), ]
rownames(phospho) = 1:nrow(phospho)

exp_prefix='ICA_phospho'
out_dir='~/documents/Segundo_Melanoma/Results/Figures'
png(filename=paste(out_dir,paste(exp_prefix,'IC_pathway.png',sep='_'),sep='/'),
    width = 1500, height = 500, units = "px", 
    pointsize = 20,
    bg = 'transparent')
ggplot(phospho,aes(x=IC,y=pathway))+
  theme_classic(base_size=20)+
  theme(axis.line=element_blank())+
  theme(axis.text.x = element_text(angle=90,vjust=0.5, size=20))+
  theme(axis.text.y = element_text(vjust=0.5, size=20))+
  geom_tile(aes(fill=logp))+
  scale_fill_gradient2(low='white',high='purple')+
  geom_text(aes(label=round(logp,digits=2)),size=6)
graphics.off()


trans = strict_ICA_GSEA_summary[strict_ICA_GSEA_summary$group=='transcriptomics' & strict_ICA_GSEA_summary$pathway %in% pathways$Transcriptomics, c('pathway', 'IC', 'logp')]
icPave = read.delim("~/Documents/Segundo_Melanoma/Results/Figures/Figure2-S2/ICA_transcriptomics_IC_Clinical_Correlation_P_Value_all.tsv", row.names = 1)
icPave = as.matrix(icPave)

icPave_plot=melt(icPave)
icPave_plot=subset(icPave_plot,
                   X1%in%names(which(apply(icPave,1,sum)>30))&X2%in%colnames(icPave)[which(apply(icPave,2,sum)>30)])
icPave_plot$value[icPave_plot$value > 100] = 100

colnames(icPave_plot)[1:2]=c('IC','clinical_feature')
icPave_plot[,1]=factor(icPave_plot[,1])
icPave_plot$IC=paste('IC',str_pad(icPave_plot$IC,3,pad='0'),sep='_')
IC = unique(icPave_plot$IC)[!unique(icPave_plot$IC) %in% trans$IC]
dum = as.data.frame(IC)
dum$pathway = NA
dum$logp = NA
dum = dum[, c("pathway", "IC", "logp")]
trans$IC = as.character(trans$IC)
trans = rbind(trans, dum)
trans = trans[order(trans$IC), ]
rownames(trans) = 1:nrow(trans)

exp_prefix='ICA_transcriptomics'
out_dir='~/documents/Segundo_Melanoma/Results/Figures'
png(filename=paste(out_dir,paste(exp_prefix,'IC_pathway.png',sep='_'),sep='/'),
    width = 1500, height = 500, units = "px", 
    pointsize = 20,
    bg = 'transparent')
ggplot(trans,aes(x=IC,y=pathway))+
  theme_classic(base_size=20)+
  theme(axis.line=element_blank())+
  theme(axis.text.x = element_text(angle=90,vjust=0.5, size=20))+
  theme(axis.text.y = element_text(vjust=0.5, size=20))+
  geom_tile(aes(fill=logp))+
  scale_fill_gradient2(low='white',high='purple')+
  geom_text(aes(label=round(logp,digits=2)),size=6)
graphics.off()






