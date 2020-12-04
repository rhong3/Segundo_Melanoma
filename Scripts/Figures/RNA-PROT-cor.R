library(readr)
# RNA-Protein correlation
proteomics <- read_csv("Data/proteomics/proteomics.csv")
proteomics = proteomics[, -c(111,112)]
transcriptomics <- read_csv("Data/transcriptomics/transcriptomics.csv")
transcriptomics = transcriptomics[,-2]
transcriptomics = transcriptomics[transcriptomics$`Gene name` %in% intersect(proteomics$`Gene name`, transcriptomics$`Gene name`), 
                                  colnames(transcriptomics) %in% intersect(colnames(proteomics), colnames(transcriptomics))]
proteomics = proteomics[proteomics$`Gene name` %in% intersect(proteomics$`Gene name`, transcriptomics$`Gene name`), 
                        colnames(proteomics) %in% intersect(colnames(proteomics), colnames(transcriptomics))]

library(dplyr)
proteomics = proteomics %>% group_by(`Gene name`) %>% mutate_each(funs(mean), ) %>% distinct
proteomics =proteomics[order(row.names(proteomics)), order(colnames(proteomics))]
proteomics = t(proteomics)
colnames(proteomics) = proteomics[1,]
proteomics = proteomics[-1,]
proteomics =proteomics[order(row.names(proteomics)), order(colnames(proteomics))]
transcriptomics =transcriptomics[order(row.names(transcriptomics)), order(colnames(transcriptomics))]
transcriptomics = t(transcriptomics)
colnames(transcriptomics) = transcriptomics[1,]
transcriptomics = transcriptomics[-1,]
transcriptomics =transcriptomics[order(row.names(transcriptomics)), order(colnames(transcriptomics))]

library(Hmisc)
corr = data.frame()
for (i in 1:ncol(proteomics)){
  m = rcorr(as.numeric(proteomics[,i]), as.numeric(transcriptomics[,i]), type='spearman')
  corr[i, 1] = colnames(proteomics)[i]
  corr[i, 2] = m$r[1,2]
  corr[i, 3] = m$P[1,2]
}
colnames(corr) = c('Gene', 'correlation coefficient', 'P-value')
corr$significant = corr$`P-value` < 0.05
corr$correlation = corr$`correlation coefficient` > 0 
pos_cor = nrow(corr[corr$correlation == TRUE,])/nrow(corr)
neg_cor = 1-pos_cor
sig_cor = nrow(corr[corr$significant == TRUE,])/nrow(corr)

corr$correlation = gsub(TRUE, '97% positive correlation', corr$correlation)
corr$correlation = gsub(FALSE, '3% negative correlation', corr$correlation)

sig = corr[corr$significant == TRUE,]
sig$correlation = '84% significant correlation (p<0.05)'

new_cor = rbind(corr, sig)

library(ggplot2)
library(ggpubr)
library(gridExtra)
# Change histogram plot fill colors by groups
p = ggplot(new_cor, aes(x=`correlation coefficient`, fill=correlation, color=correlation)) +
  geom_histogram(position="identity", alpha=0.5, binwidth=0.01)+scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(title="mRNA and protein Spearman correlation")+ 
  theme_classic()+
  theme(legend.position="top", plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.text=element_text(size=20), axis.title = element_text(size = 20), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 16)) 
pdf(file=paste("~/documents/Segundo_Melanoma/Results/mRNA-protein-corr.pdf", sep=''),
    width=12,height=8)
grid.arrange(p,nrow=1, ncol=1)
dev.off()

write_csv(corr, '~/documents/Segundo_Melanoma/Results/mRNA-protein-corr.csv')


# Protein-Phospho correlation
proteomics <- read_csv("Data/proteomics/proteomics.csv")
proteomics = proteomics[, -c(112,113)]
phospho <- read_csv("Data/phospho/phospho.csv")
phospho = phospho[,-c(1,3,4)]
phospho = phospho[phospho$Accession %in% intersect(proteomics$Accession, phospho$Accession), 
                                  colnames(phospho) %in% intersect(colnames(proteomics), colnames(phospho))]
proteomics = proteomics[proteomics$Accession %in% intersect(proteomics$Accession, phospho$Accession), 
                        colnames(proteomics) %in% intersect(colnames(proteomics), colnames(phospho))]

library(dplyr)
phospho = phospho %>% group_by(`Accession`) %>% mutate_each(funs(mean), ) %>% distinct
phospho =phospho[order(row.names(phospho)), order(colnames(phospho))]
phospho = t(phospho)
colnames(phospho) = phospho[1,]
phospho = phospho[-1,]
phospho =phospho[order(row.names(phospho)), order(colnames(phospho))]
proteomics =proteomics[order(row.names(proteomics)), order(colnames(proteomics))]
proteomics = t(proteomics)
colnames(proteomics) = proteomics[1,]
proteomics = proteomics[-1,]
proteomics =proteomics[order(row.names(proteomics)), order(colnames(proteomics))]

library(Hmisc)
corr = data.frame()
for (i in 1:ncol(proteomics)){
  m = rcorr(as.numeric(proteomics[,i]), as.numeric(phospho[,i]), type='spearman')
  corr[i, 1] = colnames(proteomics)[i]
  corr[i, 2] = m$r[1,2]
  corr[i, 3] = m$P[1,2]
}
colnames(corr) = c('Accession', 'correlation coefficient', 'P-value')
corr$significant = corr$`P-value` < 0.05
corr$correlation = corr$`correlation coefficient` > 0 
pos_cor = nrow(corr[corr$correlation == TRUE,])/nrow(corr)
neg_cor = 1-pos_cor
sig_cor = nrow(corr[corr$significant == TRUE,])/nrow(corr)

corr$correlation = gsub(TRUE, '99% positive correlation', corr$correlation)
corr$correlation = gsub(FALSE, '1% negative correlation', corr$correlation)

sig = corr[corr$significant == TRUE,]
sig$correlation = '94% significant correlation (p<0.05)'

new_cor = rbind(corr, sig)


library(ggplot2)
library(ggpubr)
library(gridExtra)
# Change histogram plot fill colors by groups
p = ggplot(new_cor, aes(x=`correlation coefficient`, fill=correlation, color=correlation)) +
  geom_histogram(position="identity", alpha=0.5, binwidth=0.01)+scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(title="protein and phosphosite Spearman correlation")+ 
  theme_classic()+
  theme(legend.position="top", plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.text=element_text(size=20), axis.title = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) 
pdf(file=paste("~/documents/Segundo_Melanoma/Results/protein-phospho-corr.pdf", sep=''),
    width=20,height=8)
grid.arrange(p,nrow=1, ncol=1)
dev.off()

write_csv(corr, '~/documents/Segundo_Melanoma/Results/protein-phospho-corr.csv')



