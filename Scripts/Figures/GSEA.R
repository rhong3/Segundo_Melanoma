## ICA-ranked GSEA
#' 
#' Created on 1/14/2020
#' 
#' @author: RH
#' 

library(BiocManager)
library('org.Hs.eg.db')
library("reactome.db")
library("ReactomePA")
library("fgsea")
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(circlize)
library("RColorBrewer")

todolist = function(ICAfile, mix_file, outfile){
  tdlist = list()
  i = 1
  ICAfile=read.table(file=ICAfile,
                      header=T,sep='\t',row.names = 1)
  mix = read.delim(mix_file)
  rownames(mix) = mix$X
  rownames(ICAfile) = rownames(mix)
  for (feature in colnames(ICAfile)){
    sig_ICs = rownames(ICAfile[ICAfile[feature] > 30, ])
    for (IC in sig_ICs){
      tdlist[[i]] = c(IC, feature)
      i = i+1
    }
  }
  expe = data.frame(t(data.frame(tdlist)))
  colnames(expe) = c('IC', 'Feature')
  expe$Feature = gsub("lymphnode..", "lymphnode.percentage", expe$Feature)
  expe$Feature = gsub("tumor..", "tumor.percentage", expe$Feature)
  expe$Feature = gsub("tumor.percentageell.size.average", "tumor.cell.size.average", expe$Feature)
  write.csv(expe, outfile, row.names = FALSE)
}


GSEA = function(centroid_file, td_list_file, outdir){
  centroid <- read.csv(centroid_file)
  td_list <- read.csv(td_list_file)
  td_list = unique(td_list$IC)
  
  for (m in td_list){
    Gene_order = centroid[order(centroid[m]),]
    ENTREZID = mapIds(org.Hs.eg.db, as.character(Gene_order$Gene.name), 'ENTREZID', 'SYMBOL')
    Gene_order.temp <- setNames(as.numeric(unlist(Gene_order[m])), unname(ENTREZID))
    my_pathways <- reactomePathways(names(Gene_order.temp))
    summary(sapply(my_pathways, length))
    fgsea_reactome <- fgsea(pathways = my_pathways, 
                            stats = Gene_order.temp,
                            minSize=15,
                            maxSize=500,
                            nperm=100000)
    fgsea_reactome <- na.omit(fgsea_reactome[order(pval), ])
    fgsea_reactome$leadingEdge = as.character(fgsea_reactome$leadingEdge)
    fgsea_reactome.sig = fgsea_reactome[fgsea_reactome$padj < 0.01,]
    write.csv(fgsea_reactome, file = paste(outdir, m, ".csv", sep=''))
    # pdf("~/Documents/Segundo_Melanoma/proteomics/ICA/0403ICA/GSEA/Most_IC8_sig.pdf", paper = 'letter')
    # plotEnrichment(my_pathways[[head(fgsea_reactome, 1)$pathway]], Gene_order.8) + labs(title=head(fgsea_reactome, 1)$pathway)
    # dev.off()
    
    topPathwaysUp <- fgsea_reactome.sig[ES > 0][head(order(padj), n=20), pathway]
    topPathwaysDown <- fgsea_reactome.sig[ES < 0][head(order(padj), n=20), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    pdf(paste(outdir, m, "_sig.pdf"), width = 15, paper = 'a4r')
    plotGseaTable(my_pathways[topPathways], Gene_order.temp, fgsea_reactome, 
                  gseaParam = 0.5)
    dev.off()
  }
}


HMP = function(centroid_file, clinical_file, data_file, td_list_file, outdir, numb, nm){
  ica = read.csv(centroid_file, row.names=1)
  clinical = read.csv(clinical_file, row.names=1)
  clinical$prim.breslow.class = gsub("NM", "NA", clinical$prim.breslow.class)
  colnames(clinical) = gsub("lymphnode..", "lymphnode.percentage", colnames(clinical))
  colnames(clinical) = gsub("tumor..", "tumor.percentage", colnames(clinical))
  colnames(clinical) = gsub("tumor.percentageell.size.average", "tumor.cell.size.average", colnames(clinical))
  ori_data = read.csv(data_file, row.names=1)
  td_list <- read.csv(td_list_file)
  for (a in 1:nrow(td_list)){
    m = toString(droplevels(td_list[a, "Feature"]))
    n = toString(droplevels(td_list[a, "IC"]))
    ica = ica[order(-ica[n]), ]
    clinical = clinical[order(clinical[m]), ]
    sorted_data_all = ori_data[match(rownames(ica), rownames(ori_data)), match(rownames(clinical), colnames(ori_data))]
    sorted_data = data.matrix(sorted_data_all[-c(11:(nrow(sorted_data_all)-10)),])
    sorted_data = data.frame(rbind(t(clinical[m]), sorted_data))
    sorted_data_all = data.frame(rbind(t(clinical[m]), data.matrix(sorted_data_all)))
    sorted_data=sorted_data[,colSums(is.na(sorted_data[1,]))==0]
    sorted_data[] <- lapply(sorted_data, function(x) {
      if(is.factor(x)) as.numeric(as.character(x)) else x
    })
    sorted_data_all[] <- lapply(sorted_data_all, function(x) {
      if(is.factor(x)) as.numeric(as.character(x)) else x
    })
    if (numb == 1){
      sorted_data_out = merge(ica[,1:2], sorted_data[2:nrow(sorted_data), ], by=0, sort=FALSE)
      sorted_data_out = sorted_data_out[,-3]
      sorted_data_all_out = merge(ica[,1:2], sorted_data_all[2:nrow(sorted_data_all),], by=0, sort=FALSE)
      sorted_data_all_out = sorted_data_all_out[,-3]
    }else{
      sorted_data_out = merge(ica[,1:numb], sorted_data[2:nrow(sorted_data), ], by=0, sort=FALSE)
      sorted_data_all_out = merge(ica[,1:numb], sorted_data_all[2:nrow(sorted_data_all),], by=0, sort=FALSE)
    }
    
    sorted_data_out = rbind.fill(sorted_data[1,], sorted_data_out)
    sorted_data_all_out = rbind.fill(sorted_data_all[1,], sorted_data_all_out)
    colnames(sorted_data_out) = gsub("Row.names", nm, colnames(sorted_data_out))
    colnames(sorted_data_all_out) = gsub("Row.names", nm, colnames(sorted_data_all_out))
    write.csv(sorted_data_out, paste(outdir, n, "_", m, "_lite.csv", sep=""), row.names=FALSE)
    write.csv(sorted_data_all_out, paste(outdir, n, "_", m, ".csv", sep=""), row.names=FALSE)
    
    sorted_data_out[nm][1,] = m
    rownames(sorted_data_out) = unname(unlist(sorted_data_out[nm]))
    
    gn = rowAnnotation(Gene.name = anno_text(sorted_data_out$Gene.name[2:length(sorted_data_out$Gene.name)],
                                             location = 0.5, just = "center"))
    sorted_data_out = sorted_data_out %>% dplyr::select(matches("MM"))
    breaksList = seq(min(ori_data), max(ori_data), by=1)
    col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
    col_fun = colorRamp2(c(min(as.numeric(sorted_data_out[1, ])), mean(as.numeric(sorted_data_out[1, ])), max(as.numeric(sorted_data_out[1, ]))), c("blue", "white", "red"))
    anno = HeatmapAnnotation(Feature = as.numeric(sorted_data_out[1, ]), col=list(Feature=col_fun), annotation_legend_param = list(direction = "horizontal"))
    pdf(paste(outdir, n, "_", m, "_HM.pdf", sep=""), height = 5, width = 20)
    hp = Heatmap(as.matrix(sorted_data_out[2:nrow(sorted_data_out), ]), col = col, column_title = paste(n, ' vs ',m), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
                 cluster_rows = FALSE, cluster_columns = FALSE, row_split = c(rep('first 10 proteins',10), rep('last 10 proteins',10)), name = "Value", heatmap_legend_param = list(direction = "horizontal"))
    draw(hp, heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom", merge_legend = TRUE,)
    dev.off()
  }
}


# proteomics
todolist("~/documents/Segundo_Melanoma/Results/proteomics/ICA/0.0005/ICA_proteomics_IC_Clinical_Correlation_P_Value_all.tsv",
         "~/documents/Segundo_Melanoma/Results/proteomics/ICA/ICA_proteomics_IC_mean_mixing_score.txt",
         "~/documents/Segundo_Melanoma/Results/proteomics/ICA/0.0005/significant_IC_clinical.csv")
GSEA("~/documents/Segundo_Melanoma/Results/proteomics/ICA/MG_ICA_proteomics_IC_centroid.csv", 
     "~/documents/Segundo_Melanoma/Results/proteomics/ICA/0.0005/significant_IC_clinical.csv", 
     "~/documents/Segundo_Melanoma/Results/proteomics/GSEA/relax/")
HMP(centroid_file = "~/documents/Segundo_Melanoma/Results/proteomics/ICA/MG_ICA_proteomics_IC_centroid.csv",
    clinical_file = "~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv",
    data_file = "~/documents/Segundo_Melanoma/Data/proteomics/ICA_proteomics.csv",
    td_list_file = "~/documents/Segundo_Melanoma/Results/proteomics/ICA/0.0005/significant_IC_clinical.csv",
    outdir = "~/documents/Segundo_Melanoma/Results/proteomics/GSEA/relax/", numb = 2, nm = "Accession")

# transcriptomics
todolist("~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/0.00001/ICA_transcriptomics_IC_Clinical_Correlation_P_Value_all.tsv",
         "~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/ICA_transcriptomics_IC_mean_mixing_score.txt",
         "~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/0.00001/significant_IC_clinical.csv")
GSEA("~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/MG_ICA_transcriptomics_IC_centroid.csv", 
     "~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/0.00001/significant_IC_clinical.csv", 
     "~/documents/Segundo_Melanoma/Results/transcriptomics/GSEA/strict/")
HMP(centroid_file = "~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/MG_ICA_transcriptomics_IC_centroid.csv",
    clinical_file = "~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv",
    data_file = "~/documents/Segundo_Melanoma/Data/transcriptomics/ICA_transcriptomics.csv",
    td_list_file = "~/documents/Segundo_Melanoma/Results/transcriptomics/ICA/0.0005/significant_IC_clinical.csv",
    outdir = "~/documents/Segundo_Melanoma/Results/transcriptomics/GSEA/relax/", numb = 1, nm="Gene.name")

# phospho
todolist("~/documents/Segundo_Melanoma/Results/phospho/ICA/0.00001/ICA_phospho_IC_Clinical_Correlation_P_Value_all.tsv",
         "~/documents/Segundo_Melanoma/Results/phospho/ICA/ICA_phospho_IC_mean_mixing_score.txt",
         "~/documents/Segundo_Melanoma/Results/phospho/ICA/0.00001/significant_IC_clinical.csv")
GSEA("~/documents/Segundo_Melanoma/Results/phospho/ICA/MG_ICA_phospho_IC_centroid.csv", 
     "~/documents/Segundo_Melanoma/Results/phospho/ICA/0.00001/significant_IC_clinical.csv", 
     "~/documents/Segundo_Melanoma/Results/phospho/GSEA/strict/")
HMP(centroid_file = "~/documents/Segundo_Melanoma/Results/phospho/ICA/MG_ICA_phospho_IC_centroid.csv",
    clinical_file = "~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv",
    data_file = "~/documents/Segundo_Melanoma/Data/phospho/ICA_phospho.csv",
    td_list_file = "~/documents/Segundo_Melanoma/Results/phospho/ICA/0.0005/significant_IC_clinical.csv",
    outdir = "~/documents/Segundo_Melanoma/Results/phospho/GSEA/relax/", numb = 3, nm="Modified_sequence")



# # Proteomics subtypes
# todolist("~/documents/Segundo_Melanoma/Results/proteomics/ICA/subtype/0.0005/ICA_proteomics_IC_Clinical_Correlation_P_Value_all.tsv",
#          "~/documents/Segundo_Melanoma/Results/proteomics/ICA/ICA_proteomics_IC_mean_mixing_score.txt",
#          "~/documents/Segundo_Melanoma/Results/proteomics/ICA/subtype/0.0005/significant_IC_clinical.csv")
# 

