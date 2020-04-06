# Make density plots
library(ggplot2)
library(reshape2)
library(dplyr)
prot <- read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
phos <- read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")

# prot <- read.csv("/beegfs/rh2740/Fenyo/Lund_Melanoma/Data/proteomics/proteomics.csv")
# phos <- read.csv("/beegfs/rh2740/Fenyo/Lund_Melanoma/Data/phospho/phospho.csv")
rownames(prot) = prot$Accession
rownames(phos) = phos$Modified_sequence
prot = prot[,-c(112:113)]
phos = phos[,-c(2:4)]

prot.melt = melt(prot, id.vars = "Accession")
p <- ggplot(aes(x=value, colour=variable), data=prot.melt)
p + geom_density()+theme(panel.background = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), plot.background = element_blank())+xlim(-8, 8)
ggsave('~/documents/Segundo_Melanoma/Results/prot_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)

phos.melt = melt(phos, id.vars = "Modified_sequence")
h <- ggplot(aes(x=value, colour=variable), data=phos.melt)
h + geom_density() + theme(panel.background = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), plot.background = element_blank())+xlim(-8, 8)
ggsave('~/documents/Segundo_Melanoma/Results/phos_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)

phosx = read.delim("~/Documents/Segundo_Melanoma/Data/phospho/MM_Phosho_DIA_lund_Cohort_total_peptides_NORMALIZED.txt")
phosx = phosx[,-c(2:8)]
colnames(phosx)[1] = 'Modified_sequence' 
phosx$MM790 = rowMeans(phosx[c('MM790.1.', 'MM790.2.')], na.rm=TRUE)
phosx$MM807 = rowMeans(phosx[c('MM807.1.', 'MM807.2.')], na.rm=TRUE)
phosx$MM808 = phosx$MM808_LG
phosx = select(phosx, -c(MM790.1., MM790.2., MM807.1., MM807.2., MM808_LG))
Phospho_his_clin <- read.csv("~/Documents/Segundo_Melanoma/Data/clinical&histology.csv", row.names=1)
phosx.sample = subset(phosx , select = names(phosx) %in% c(intersect(colnames(phosx), rownames(Phospho_his_clin)), "Modified_sequence", "Accession", "Description", "Gene name"))
phosx.sample = na.omit(phosx.sample)
rownames(phosx.sample) = phosx.sample$Modified_sequence
phosx.melt = melt(phosx.sample, id.vars = "Modified_sequence")
h <- ggplot(aes(x=value, colour=variable), data=phosx.melt)
h + geom_density() + theme(panel.background = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), plot.background = element_blank())+xlim(-8, 8)
ggsave('~/documents/Segundo_Melanoma/Results/RAW_phos_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)
