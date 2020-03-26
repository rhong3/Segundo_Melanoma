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
    colnames(summ1) = c("from", "to", "-logp")
    colnames(summ2) = c("from", "to", "-logp")
    
    xxx = xxxx[xxxx$V2.y %in% unqp$name & xxxx$V2.x %in% unqp$name,]
    xxx = xxx[, c(2,1)]
    colnames(xxx) = c('from', 'to')
    xxx['-logp'] = mean(summ1$`-logp`)
    
    wxxx = xxx[ ,c(1,2)]
    colnames(wxxx) = c('from', 'to2')
    summ2 = left_join(summ2, wxxx, by='from')
    summ2[, 4] = as.character(summ2[, 4])
    for (row in 1:nrow(summ2)){
      if (is.na(summ2[row, 4])){
        summ2[row, 4] = summ2[row, 1]
      }
    }
    
    summ2 = summ2[, c(4,2,3)]
    colnames(summ2) = c('from', 'to', '-logp')
    unqp = distinct(summ2['from'])
    colnames(unqp) = c('name')
    unqp['group'] = 'pathway'
    unq = rbind(unqp, unqc, unqi)
    unq = unique(unq)
    
    summall = rbind(summ1, summ2, xxx)
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
    colrs <- c("tomato", "gold", "light blue")
    V(g)$group = gsub("pathway", "light blue", V(g)$group)
    V(g)$group = gsub("clinical", "tomato", V(g)$group)
    V(g)$group = gsub("IC", "gold", V(g)$group)
    V(g)$color <- V(g)$group
    deg <- degree(g, mode="all")
    V(g)$size <- 2*log(deg)
    E(g)$width <- E(g)$`-logp`/3
    E(g)$arrow.size <- .5
    E(g)$edge.color <- "gray80"
    # V(g)$label <- NA
    edge.start <- get.edges(g, 1:ecount(g))[,2] 
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

# For multiomics validated pathways
# Pathway, IC, clinical
Relation <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("~/Documents/Segundo_Melanoma/Results/ReactomePathways.txt", header=FALSE)
Pathways = Pathways[Pathways['V3'] == 'Homo sapiens', c(1,2)]
xxxx = na.omit(left_join(Relation, Pathways, by="V1"))[, c(2,3)]
colnames(xxxx) = c('V1', 'V2')
xxxx = na.omit(left_join(xxxx, Pathways, by="V1"))[, c(2,3)]

for (a in c('strict', 'relax')){
  summ = read.csv(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary_joined.csv", sep=''))
  summ$IC.prot = gsub('IC', 'prot', summ$IC.prot)
  summ$IC.trans = gsub('IC', 'trans', summ$IC.trans)
  summ$IC.phos = gsub('IC', 'phos', summ$IC.phos)
  summ['-logp.prot'] = -log(summ$padj.prot)
  summ['-logp.trans'] = -log(summ$padj.trans)
  summ['-logp.phos'] = -log(summ$padj.phos)
  summ = summ[c('pathway', 'clinical', 'IC.prot', 'IC.trans', 'IC.phos', '-logp.prot', '-logp.trans', '-logp.phos')]
  IC2clinp = summ[c('IC.prot', 'clinical', '-logp.prot')]
  IC2clint = summ[c('IC.trans', 'clinical', '-logp.trans')]
  IC2clinh = summ[c('IC.phos', 'clinical', '-logp.phos')]
  colnames(IC2clinp) = c('from', 'to', '-logp')
  colnames(IC2clint) = c('from', 'to', '-logp')
  colnames(IC2clinh) = c('from', 'to', '-logp')
  IC2clin = rbind(IC2clinp, IC2clint, IC2clinh)
  IC2clin = distinct(IC2clin)
  
  pw2ICp = summ[c('pathway', 'IC.prot', '-logp.prot')]
  pw2ICt = summ[c('pathway', 'IC.trans', '-logp.trans')]
  pw2ICh = summ[c('pathway', 'IC.phos', '-logp.phos')]
  colnames(pw2ICp) = c('from', 'to', '-logp')
  colnames(pw2ICt) = c('from', 'to', '-logp')
  colnames(pw2ICh) = c('from', 'to', '-logp')
  pw2IC = rbind(pw2ICp, pw2ICt, pw2ICh)
  pw2IC = distinct(pw2IC)
  
  unqp = unique(pw2IC['from'])
  xxx = xxxx[xxxx$V2.y %in% unqp$from & xxxx$V2.x %in% unqp$from,]
  xxx = xxx[, c(2,1)]
  colnames(xxx) = c('from', 'to')
  xxx['-logp'] = mean(pw2IC$`-logp`)
  wxxx = xxx[, c(1,2)]
  colnames(wxxx) = c('from', 'to2')
  pw2ICM = left_join(pw2IC, wxxx, by='from')
  
  pw2ICM[, 4] = as.character(pw2ICM[, 4])
  for (row in 1:nrow(pw2ICM)){
    if (is.na(pw2ICM[row, 4])){
      pw2ICM[row, 4] = pw2ICM[row, 1]
    }
  }
  pw2IC = pw2ICM[, c(4,2,3)]
  colnames(pw2IC) = c('from', 'to', '-logp')
  pw2IC = distinct(pw2IC)
  
  ver = distinct(rbind(IC2clin, pw2IC, xxx))
  unqcl = unique(IC2clin['to'])
  unqcl['group'] = 'clinical'
  colnames(unqcl) = c('from', 'group')
  unqIC = unique(IC2clin['from'])
  unqIC['group'] = 'IC'
  unqp1 = unique(pw2IC['from'])
  unqp2 = unique(xxx['from'])
  unqp = rbind(unqp1, unqp2)
  unqp['group'] = 'pathway'
  unq = rbind(unqcl, unqIC, unqp)
  colnames(unq) = c('name', 'group')
  unq = distinct(unq)
  
  g <- graph_from_data_frame(ver, directed=TRUE, vertices=unq)
  g <- simplify(g, remove.multiple = F, remove.loops = T, edge.attr.comb="sum")
  colrs <- c("tomato", "gold", "light blue")
  V(g)$group = gsub("pathway", "light blue", V(g)$group)
  V(g)$group = gsub("clinical", "tomato", V(g)$group)
  V(g)$group = gsub("IC", "gold", V(g)$group)
  V(g)$color <- V(g)$group
  deg <- degree(g, mode="all")
  V(g)$size <- 2*log(deg)
  E(g)$width <- E(g)$`-logp`/3
  E(g)$arrow.size <- .5
  E(g)$edge.color <- "gray80"
  # V(g)$label <- NA
  edge.start <- get.edges(g, 1:ecount(g))[,2] 
  edge.col <- V(g)$color[edge.start]
  l = layout.kamada.kawai(g)
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  
  pdf(paste("~/Documents/Segundo_Melanoma/Results/", a, "_ICA_GSEA_summary_joined.pdf", sep=''), height = 20, width = 20)
  plot(g, edge.color=edge.col, edge.curved=.2, vertex.label.color="black", rescale=F, layout=l*1, 
       vertex.label.cex=0.8, main=paste(a, ' joined ICA GSEA', sep=''), vertex.label.degree=pi/2)
  legend(x=-1, y=-.8, c("clinical features", "independent components", "pathways"), pch=21,
         col="#777777", pt.bg=colrs, pt.cex=2, cex=2, bty="n", ncol=1)
  dev.off()
}  
