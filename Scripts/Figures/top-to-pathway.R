# 10 genes pathway involvement
toplist = c('ADAM10', 'SCAI', 'TEX30', 'HMOX1', 'CDK4', 'CTNND1', 'DDX11', 'FGA', 'PAEP', 'PIK3CB')

library("org.Hs.eg.db")
library("reactome.db")
annotation <- select(org.Hs.eg.db, keys=toplist, columns=c('SYMBOL', 'ENTREZID'), keytype="SYMBOL")
annotation.path <- select(reactome.db, keys = unique(annotation$ENTREZID), columns=c("PATHID", "PATHNAME"), keytype='ENTREZID')
anno = merge(annotation, annotation.path, by='ENTREZID')
anno$PATHNAME = gsub('Homo sapiens: ', '', anno$PATHNAME)

no_pathway = c('TEX30', 'PAEP')

relax <- read_csv("Results/relax_ICA_GSEA_summary.csv")
strict <- read_csv("Results/strict_ICA_GSEA_summary.csv")

for (m in unique(anno$SYMBOL)){
  temp = anno[anno$SYMBOL == m,]
  relax[m] = as.numeric(relax$pathway %in% temp$PATHNAME)
  strict[m] = as.numeric(strict$pathway %in% temp$PATHNAME)
}


relax['toplists_num'] = rowSums(relax[,12:21])
strict['toplists_num'] = rowSums(strict[,12:21])
relax = relax[order(relax$toplists_num, decreasing = TRUE),]
strict = strict[order(strict$toplists_num, decreasing = TRUE),]
relax = relax[relax$toplists_num > 0, ]
strict = strict[strict$toplists_num > 0, ]

for (i in 1:nrow(relax)){
  temp = read_csv(paste("Results/", relax[i, 'group'], '/GSEA/relax/', relax[i, 'IC'], '_', relax[i, 'clinical'], '_lite.csv', sep=''))
  inter = intersect(temp$Gene.name, toplist)
  for (f in inter){
    relax[i, f] = relax[i, f]*100
  }
}
relax$importance_score = rowSums(relax[,12:21])


for (i in 1:nrow(strict)){
  temp = read_csv(paste("Results/", strict[i, 'group'], '/GSEA/strict/', strict[i, 'IC'], '_', strict[i, 'clinical'], '_lite.csv', sep=''))
  inter = intersect(temp$Gene.name, toplist)
  for (f in inter){
    strict[i, f] = strict[i, f]*100
  }
}
strict$importance_score = rowSums(strict[,12:21])

relax = relax[order(relax$importance_score, decreasing = TRUE),]
strict = strict[order(strict$importance_score, decreasing = TRUE),]

write.csv(relax, "Results/toplist_relax_ICA_GSEA_summary.csv", row.names = FALSE)
write.csv(strict, "Results/toplist_strict_ICA_GSEA_summary.csv", row.names = FALSE)




