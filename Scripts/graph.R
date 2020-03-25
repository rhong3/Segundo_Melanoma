# Pathway iGraph
library(igraph)
library(dplyr)

# Pathway, IC, clinical
Relation <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathways.txt", header=FALSE)
Pathways = Pathways[Pathways['V3'] == 'Homo sapiens', c(1,2)]
xxxx = na.omit(left_join(Relation, Pathways, by="V1"))[, c(2,3)]
colnames(xxxx) = c('V1', 'V2')
xxxx = na.omit(left_join(xxxx, Pathways, by="V1"))[, c(2,3)]

for (a in c('strict', 'relax')){
  summ = read.csv(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary.csv", sep=''))
  summ['-logp'] = -log(summ$padj)
  for (b in c('proteomics', 'transcriptomics', 'phospho')){
    summ1 = summ[summ$group == b, c(10,11,12)]
    summ2 = summ[summ$group == b, c(1,10,12)]
    
    unqi = distinct(summ1['IC'])
    colnames(unqi) = c('name')
    unqi['group'] = 'IC'
    unqp = distinct(summ2['pathway'])
    colnames(unqp) = c('name')
    unqp['group'] = 'pathway'
    unqc = distinct(summ1['clinical'])
    colnames(unqc) = c('name')
    unqc['group'] = 'clinical'
    unq = rbind(unqp, unqc, unqi)
    colnames(summ1) = c("from", "to", "-logp")
    colnames(summ2) = c("from", "to", "-logp")
    
    xxx = xxxx[xxxx$V2.y %in% unqp$name & xxxx$V2.x %in% unqp$name,]
    xxx = xxx[, c(2,1)]
    colnames(xxx) = c('from', 'to')
    xxx['-logp'] = mean(summ1$`-logp`)
    
    summall = rbind(summ1, summ2, xxx)
    summall = unique(summall)
    g <- graph_from_data_frame(summall, directed=TRUE, vertices=unq)
    g <- simplify(g, remove.multiple = F, remove.loops = T)
    colrs <- c("tomato", "gold", "light blue")
    V(g)$group = gsub("pathway", "light blue", V(g)$group)
    V(g)$group = gsub("clinical", "tomato", V(g)$group)
    V(g)$group = gsub("IC", "gold", V(g)$group)
    V(g)$color <- V(g)$group
    deg <- degree(g, mode="all")
    V(g)$size <- 2*log(deg)
    E(g)$width <- E(g)$`-logp`/3
    E(g)$arrow.size <- .05
    E(g)$edge.color <- "gray80"
    # V(g)$label <- NA
    edge.start <- get.edges(g, 1:ecount(g))[,1] 
    edge.col <- V(g)$color[edge.start]
    l = layout.kamada.kawai(g)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    
    pdf(paste("~/Documents/Segundo_Melanoma/Results/", b, '_', a, "_ICA_GSEA_summary.pdf", sep=''), height = 20, width = 20)
    plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1, 
         vertex.label.cex=0.6, main=paste(b, '(', a, ') ICA GSEA', sep=''), vertex.label.degree=pi/2)
    legend(x=-1, y=-.8, c("clinical features", "independent components", "pathways"), pch=21,
           col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
    dev.off()
  }
}
