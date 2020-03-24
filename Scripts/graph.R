# Pathway iGraph
library(igraph)
library(dplyr)

strict = read.csv("~/Documents/Segundo_Melanoma/Results/strict_ICA_GSEA_summary.csv")
relax = read.csv("~/Documents/Segundo_Melanoma/Results/relax_ICA_GSEA_summary.csv")
strict['-logp'] = -log(strict$padj)
relax['-logp'] = -log(relax$padj)
strict.prot = strict[strict$group == 'proteomics', c(1,7,11,12)]
unqp = distinct(strict.prot['pathway'])
colnames(unqp) = c('name')
unqp['group'] = 'pathway'
unqc = distinct(strict.prot['clinical'])
colnames(unqc) = c('name')
unqc['group'] = 'clinical'
unq = rbind(unqp, unqc)
colnames(strict.prot) = c("from", "size", "to", "-logp")
strict.prot = strict.prot[,c(1,3,2,4)]
g <- graph_from_data_frame(strict.prot, directed=TRUE, vertices=unq)
g <- simplify(g, remove.multiple = F, remove.loops = T)
