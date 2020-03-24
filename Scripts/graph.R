# Pathway iGraph
library(igraph)
library(dplyr)

strict = read.csv("~/Documents/Segundo_Melanoma/Results/strict_ICA_GSEA_summary.csv")
relax = read.csv("~/Documents/Segundo_Melanoma/Results/relax_ICA_GSEA_summary.csv")

strict['-logp'] = -log(strict$padj)
relax['-logp'] = -log(relax$padj)
strict.prot1 = strict[strict$group == 'proteomics', c(10,11,12)]
strict.prot2 = strict[strict$group == 'proteomics', c(1,10,12)]

unqi = distinct(strict.prot1['IC'])
colnames(unqi) = c('name')
unqi['group'] = 'IC'
unqp = distinct(strict.prot2['pathway'])
colnames(unqp) = c('name')
unqp['group'] = 'pathway'
unqc = distinct(strict.prot1['clinical'])
colnames(unqc) = c('name')
unqc['group'] = 'clinical'
unq = rbind(unqp, unqc, unqi)
colnames(strict.prot1) = c("from", "to", "-logp")
colnames(strict.prot2) = c("from", "to", "-logp")
strict.prot = rbind(strict.prot1, strict.prot2)
strict.prot = unique(strict.prot)
g <- graph_from_data_frame(strict.prot, directed=TRUE, vertices=unq)
g <- simplify(g, remove.multiple = F, remove.loops = F)
colrs <- c("tomato", "gold", "light blue")
V(g)$group = gsub("pathway", "light blue", V(g)$group)
V(g)$group = gsub("clinical", "tomato", V(g)$group)
V(g)$group = gsub("IC", "gold", V(g)$group)
V(g)$color <- V(g)$group
deg <- degree(g, mode="all")
V(g)$size <- 2*log(deg)
E(g)$width <- E(g)$`-logp`/6
E(g)$arrow.size <- .05
E(g)$edge.color <- "gray80"
# V(g)$label <- NA
edge.start <- get.edges(g, 1:ecount(g))[,1] 
edge.col <- V(g)$color[edge.start]
l = layout.kamada.kawai(g)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1.5, vertex.label.cex=.35)
legend(x=-1.8, y=-1.3, c("clinical features", "independent components", "pathways"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

