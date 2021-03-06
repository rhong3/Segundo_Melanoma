# Pathway iGraph lite
library(igraph)
library(dplyr)

inlist = c('tumor.cell.size.average', 'conncetive.tissue.average', 'predominant.tumor.cell.shape.average',
              'Average.necrosis', 'Average.tumor.percentage', 'survival_3yr',
              'pigment.score.average', 'predominant.cytoplasm.average', 'Average.lymphatic.score', 'prim.breslow',
           'survival_6mo', 'Local.Visceral', 'stage', 'Local.Cutaneous', 'Targeted.BRAF.treatment..vemurafenib.')
# Pathway, IC, clinical
Relation <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathways.txt", header=FALSE)
Pathways = Pathways[Pathways['V3'] == 'Homo sapiens', c(1,2)]
xxxx = na.omit(left_join(Relation, Pathways, by="V1"))[, c(2,3)]
colnames(xxxx) = c('V1', 'V2')
xxxx = na.omit(left_join(xxxx, Pathways, by="V1"))[, c(2,3)]

# For figure
a = 'relax'
b = "transcriptomics"
summm = read.csv(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary.csv", sep=''))
summm['-logp'] = -log(summm$padj)
summ = summm[summm$group == b, c(1,11,12)]

colnames(summ) = c("from", "to", "-logp")
summall = unique(summ)
summall = summall[summall$to %in% inlist,]
summall_g = summall[summall$`-logp`> 5,]

unqa = distinct(summall_g['from'])
colnames(unqa) = c('name')
unqa$group = 'pathway'
unqb = distinct(summall_g['to'])
colnames(unqb) = c('name')
unqb$group = 'clinical'
unq = rbind(unqa, unqb)




g <- graph_from_data_frame(summall_g, directed=TRUE, vertices=unq)
g <- simplify(g, remove.multiple = F, remove.loops = T, edge.attr.comb="sum")
colrs <- c("tomato", "light blue")
V(g)$group = gsub("pathway", "light blue", V(g)$group)
V(g)$group = gsub("clinical", "tomato", V(g)$group)
V(g)$color <- V(g)$group

V(g)$label = ifelse(nchar(V(g)$name) < 12| V(g)$group == "tomato", V(g)$name, NA)


deg <- degree(g, mode="all")
V(g)$size <- sqrt(deg)
E(g)$width <- E(g)$`-logp`/5
E(g)$arrow.size <- .3
E(g)$edge.color <- "gray80"
# V(g)$label <- NA
edge.start <- get.edges(g, 1:ecount(g))[,2]
edge.col <- V(g)$color[edge.start]
l = layout_with_kk(g)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf(paste("~/Documents/Segundo_Melanoma/Results/graph/Figure/", b, '_', a, "_ICA_GSEA_summary.pdf", sep=''), height = 20, width = 22)
plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=T, layout=l*1,
     vertex.label.cex=2.5, main=paste(b, '(', a, ') ICA GSEA', sep=''), vertex.label.degree=pi/4, vertex.label.dist=0)
legend(x=-1, y=-0.8, c("features", "pathways"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
dev.off()




# Old
for (a in c('strict', 'relax')){
  summm = read.csv(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary.csv", sep=''))
  # summm = summm[(summm$clinical %in% histology), ]
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
    if (a == 'relax'){
      if (b == 'phospho'){ma = 18}
      else if(b == 'proteomics'){ma = 15}
      else{ma = 23}
      V(g)$label = ifelse(degree(g) > ma | V(g)$group == "tomato", V(g)$name, NA)
    }else{
      if (b == 'phospho'){ma = 18}
      else if(b == 'proteomics'){ma = 10}
      else{ma = 3}
      V(g)$label = ifelse(degree(g) > ma | V(g)$group == "tomato", V(g)$name, NA)
    }
    
    deg <- degree(g, mode="all")
    V(g)$size <- log(deg)
    E(g)$width <- E(g)$`-logp`/5
    E(g)$arrow.size <- .3
    E(g)$edge.color <- "gray80"
    # V(g)$label <- NA
    edge.start <- get.edges(g, 1:ecount(g))[,2] 
    edge.col <- V(g)$color[edge.start]
    l = layout_in_circle(g)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    
    pdf(paste("~/Documents/Segundo_Melanoma/Results/graph/No_IC/", b, '_', a, "_ICA_GSEA_summary.pdf", sep=''), height = 20, width = 22)
    plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1,
         vertex.label.cex=1.25, main=paste(b, '(', a, ') ICA GSEA', sep=''), vertex.label.degree=pi/4, vertex.label.dist=0)
    legend(x=-1, y=-0.8, c("features", "pathways"), pch=21,
           col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
    dev.off()
  }
}


# For multiomics validated pathways
# Pathway, IC, clinical
histology = c('tumor.cell.size.average', 'conncetive.tissue.average', 'predominant.tumor.cell.shape.average',
              'Average.necrosis', 'lymphocyte.density.average', 'Average.necrosis', 'Average.tumor.percentage',
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
  # summ = summ[(summ$clinical %in% histology), ]
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
  V(g)$label = V(g)$name
  V(g)$color <- V(g)$group 
  deg <- degree(g, mode="all")
  V(g)$size <- log(deg)
  E(g)$width <- E(g)$`-logp`/5
  E(g)$arrow.size <- .3
  E(g)$edge.color <- "gray80"
  edge.start <- get.edges(g, 1:ecount(g))[,2] 
  edge.col <- V(g)$color[edge.start]
  l = layout_in_circle(g)
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  
  pdf(paste("~/Documents/Segundo_Melanoma/Results/graph/No_IC/", a, "_ICA_GSEA_summary_joined.pdf", sep=''), height = 20, width = 22)
  plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1, vertex.label = ifelse(degree(g) > 5 | V(g)$group == "tomato", V(g)$label, NA),
       vertex.label.cex=1.25, main=paste(a, ' joined ICA GSEA', sep=''), vertex.label.degree=pi/4, vertex.label.dist=0)
  legend(x=-1, y=-0.8, c("features", "pathways"), pch=21,
         col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
  dev.off()
}  
