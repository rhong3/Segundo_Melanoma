# Pathway iGraph lite
library(igraph)
library(dplyr)

histology = c('tumor.cell.size.average', 'conncetive.tissue.average', 'predominant.tumor.cell.shape.average',
              'Average.necrosis', 'lymphocyte.density.average', 'Average.necrosis',
              'pigment.score.average', 'predominant.cytoplasm.average', 'lymphocyte.distribution.average',
              'Average.lymphatic.score', 'average.adjacent.lymphnode.percentage')
# Pathway, IC, clinical
Relation <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathways.txt", header=FALSE)
Pathways = Pathways[Pathways['V3'] == 'Homo sapiens', c(1,2)]
xxxx = na.omit(left_join(Relation, Pathways, by="V1"))[, c(2,3)]
colnames(xxxx) = c('V1', 'V2')
xxxx = na.omit(left_join(xxxx, Pathways, by="V1"))[, c(2,3)]

for (a in c('strict', 'relax')){
  summm = read.csv(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary.csv", sep=''))
  summm = summm[!(summm$clinical %in% histology), ]
  summm['-logp'] = -log(summm$padj)
  for (b in c('proteomics', 'transcriptomics', 'phospho')){
    summ = summm[summm$group == b, c(1,11,12)]
    
    unqp = distinct(summ['pathway'])
    colnames(unqp) = c('name')
    unqp['group'] = 'pathway'
    unqc = distinct(summ['clinical'])
    colnames(unqc) = c('name')
    unqc['group'] = 'clinical'
    colnames(summ) = c("from", "to", "-logp")

    
    xxx = xxxx[xxxx$V2.y %in% unqp$name & xxxx$V2.x %in% unqp$name,]
    xxx = xxx[, c(2,1)]
    colnames(xxx) = c('from', 'to')
    xxx['-logp'] = mean(summ$`-logp`)
    
    wxxx = xxx[ ,c(1,2)]
    colnames(wxxx) = c('from', 'to2')
    summ = left_join(summ, wxxx, by='from')
    summ[, 4] = as.character(summ[, 4])
    for (row in 1:nrow(summ)){
      if (is.na(summ[row, 4])){
        summ[row, 4] = summ[row, 1]
      }
    }
    
    summ = summ[, c(4,2,3)]
    colnames(summ) = c('from', 'to', '-logp')
    unqp = distinct(summ['from'])
    colnames(unqp) = c('name')
    unqp['group'] = 'pathway'
    unq = rbind(unqp, unqc)
    unq = unique(unq)
    
    summall = rbind(summ, xxx)
    summall = unique(summall)
    
    unqa = distinct(summall['from'])
    colnames(unqa) = c('name')
    unqb = distinct(summall['to'])
    colnames(unqb) = c('name')
    unqx = rbind(unqa, unqb)
    unqx= unique(unqx)
    unq = left_join(unqx, unq, by = 'name')
    unq[is.na(unq$group), 2] = 'pathway'
    
    g <- graph_from_data_frame(summall, directed=TRUE, vertices=unq)
    g <- simplify(g, remove.multiple = F, remove.loops = T, edge.attr.comb="sum")
    colrs <- c("tomato", "light blue")
    V(g)$group = gsub("pathway", "light blue", V(g)$group)
    V(g)$group = gsub("clinical", "tomato", V(g)$group)
    V(g)$color <- V(g)$group
    deg <- degree(g, mode="all")
    V(g)$size <- log(deg)
    E(g)$width <- E(g)$`-logp`/5
    E(g)$arrow.size <- .3
    E(g)$edge.color <- "gray80"
    # V(g)$label <- NA
    edge.start <- get.edges(g, 1:ecount(g))[,2] 
    edge.col <- V(g)$color[edge.start]
    l = layout.kamada.kawai(g)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    
    pdf(paste("~/Documents/Segundo_Melanoma/Results/graph/clin_", b, '_', a, "_ICA_GSEA_summary.pdf", sep=''), height = 20, width = 20)
    plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1, 
         vertex.label.cex=0.6, main=paste(b, '(', a, ') ICA GSEA', sep=''), vertex.label.degree=pi/2)
    legend(x=-1, y=-.8, c("clinical features", "pathways"), pch=21,
           col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
    dev.off()
  }
}


# For multiomics validated pathways
# Pathway, IC, clinical
histology = c('tumor.cell.size.average', 'conncetive.tissue.average', 'predominant.tumor.cell.shape.average',
              'Average.necrosis', 'lymphocyte.density.average', 'Average.necrosis',
              'pigment.score.average', 'predominant.cytoplasm.average', 'lymphocyte.distribution.average',
              'Average.lymphatic.score', 'average.adjacent.lymphnode.percentage')
# Pathway, IC, clinical
Relation <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathways.txt", header=FALSE)
Pathways = Pathways[Pathways['V3'] == 'Homo sapiens', c(1,2)]
xxxx = na.omit(left_join(Relation, Pathways, by="V1"))[, c(2,3)]
colnames(xxxx) = c('V1', 'V2')
xxxx = na.omit(left_join(xxxx, Pathways, by="V1"))[, c(2,3)]

for (a in c('strict', 'relax')){
  summ = read.csv(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary_joined.csv", sep=''))
  summ = summ[(summ$clinical %in% histology), ]
  summ['-logp'] = -(log(summ$padj.prot)+log(summ$padj.trans)+log(summ$padj.phos))/3
  summ = summ[c('pathway', 'clinical', '-logp')]
  unqp = distinct(summ['pathway'])
  colnames(unqp) = c('name')
  unqp['group'] = 'pathway'
  unqc = distinct(summ['clinical'])
  colnames(unqc) = c('name')
  unqc['group'] = 'clinical'
  colnames(summ) = c("from", "to", "-logp")
  
  
  xxx = xxxx[xxxx$V2.y %in% unqp$name & xxxx$V2.x %in% unqp$name,]
  xxx = xxx[, c(2,1)]
  colnames(xxx) = c('from', 'to')
  xxx['-logp'] = mean(summ$`-logp`)
  
  wxxx = xxx[ ,c(1,2)]
  colnames(wxxx) = c('from', 'to2')
  summ = left_join(summ, wxxx, by='from')
  summ[, 4] = as.character(summ[, 4])
  for (row in 1:nrow(summ)){
    if (is.na(summ[row, 4])){
      summ[row, 4] = summ[row, 1]
    }
  }
  
  summ = summ[, c(4,2,3)]
  colnames(summ) = c('from', 'to', '-logp')
  unqp = distinct(summ['from'])
  colnames(unqp) = c('name')
  unqp['group'] = 'pathway'
  unq = rbind(unqp, unqc)
  unq = unique(unq)
  
  summall = rbind(summ, xxx)
  summall = unique(summall)
  
  unqa = distinct(summall['from'])
  colnames(unqa) = c('name')
  unqb = distinct(summall['to'])
  colnames(unqb) = c('name')
  unqx = rbind(unqa, unqb)
  unqx= unique(unqx)
  unq = left_join(unqx, unq, by = 'name')
  unq[is.na(unq$group), 2] = 'pathway'
  
  g <- graph_from_data_frame(summall, directed=TRUE, vertices=unq)
  g <- simplify(g, remove.multiple = F, remove.loops = T, edge.attr.comb="sum")
  colrs <- c("tomato", "light blue")
  V(g)$group = gsub("pathway", "light blue", V(g)$group)
  V(g)$group = gsub("clinical", "tomato", V(g)$group)
  V(g)$color <- V(g)$group
  deg <- degree(g, mode="all")
  V(g)$size <- log(deg)
  E(g)$width <- E(g)$`-logp`/5
  E(g)$arrow.size <- .3
  E(g)$edge.color <- "gray80"
  # V(g)$label <- NA
  edge.start <- get.edges(g, 1:ecount(g))[,2] 
  edge.col <- V(g)$color[edge.start]
  l = layout.kamada.kawai(g)
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  
  pdf(paste("~/Documents/Segundo_Melanoma/Results/graph/hist_", a, "_ICA_GSEA_summary_joined.pdf", sep=''), height = 20, width = 20)
  plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1, 
       vertex.label.cex=0.8, main=paste(a, ' joined ICA GSEA', sep=''), vertex.label.degree=pi/2)
  legend(x=-1, y=-.8, c("clinical features", "pathways"), pch=21,
         col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
  dev.off()
}  
