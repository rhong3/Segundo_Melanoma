# Make density plots
library(ggplot2)
library(reshape2)
# prot <- read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
# phos <- read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")

prot <- read.csv("/beegfs/rh2740/Fenyo/Lund_Melanoma/Data/proteomics/proteomics.csv")
phos <- read.csv("/beegfs/rh2740/Fenyo/Lund_Melanoma/Data/phospho/phospho.csv")
rownames(prot) = prot$Accession
rownames(phos) = phos$Modified_sequence
prot = prot[,-c(112:113)]
phos = phos[,-c(2:4)]

prot.melt = melt(prot, id.vars = "Accession")
prot.melt = prot.melt[,c(1,3)]
p <- ggplot(aes(x=value, colour=variable), data=prot.melt)
p + geom_density()
ggsave('/beegfs/rh2740/Fenyo/Lund_Melanoma/Results/prot_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)

phos.melt = melt(phos, id.vars = "Modified_sequence")
h <- ggplot(aes(x=value, colour=variable), data=phos.melt)
h + geom_desity()
ggsave('/beegfs/rh2740/Fenyo/Lund_Melanoma/Results/phos_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)

