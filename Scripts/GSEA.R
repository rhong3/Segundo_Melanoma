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
library("pheatmap")
library(RColorBrewer)
library(dplyr)

todolist = function(ICAfile, mix_file, outfile){
  tdlist = list()
  i = 1
  ICAfile=read.table(file=ICAfile,
                      header=T,sep='\t',row.names = 1)
  mix = read.delim(mix_file)
  rownames(mix) = mix$X
  rownames(ICAfile) = rownames(mix)
  for (feature in colnames(ICAfile)){
    sig_ICs = rownames(ICAfile[ICAfile[feature] > 50, ])
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
    # pdf("~/Documents/Lund_Melanoma/proteomics/ICA/0403ICA/GSEA/Most_IC8_sig.pdf", paper = 'letter')
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
    ica = ica[order(ica[n]), ]
    clinical = clinical[order(clinical[m]), ]
    sorted_data_all = ori_data[match(rownames(ica), rownames(ori_data)), match(rownames(clinical), colnames(ori_data))]
    sorted_data = data.matrix(sorted_data_all[-c(26:(nrow(sorted_data_all)-26)),])
    sorted_data = rbind(data.frame(t(clinical[m])), sorted_data)
    sorted_data_all = rbind(data.frame(t(clinical[m])), sorted_data_all)
    fw = data.frame(as.factor(sorted_data[1, ]))
    colnames(fw) = m
    pdf(paste(outdir, n, "_", m, "_HM.pdf", sep=""), width = 20, height = 10)
    pheatmap(sorted_data[2:nrow(sorted_data), ], cluster_cols = F, cluster_rows = F, annotation_col = fw, legend = FALSE, 
             annotation_legend = TRUE, fontsize_col = 6, main = paste(n, ' vs ',m), show_rownames = T)
    dev.off()
    if (numb == 1){
      sorted_data_out = merge(ica[,1:2], sorted_data[2:nrow(sorted_data), ], by=0)
      sorted_data_out = sorted_data_out[,-3]
      sorted_data_all_out = merge(ica[,1:2], sorted_data_all[2:nrow(sorted_data_all),], by=0)
      sorted_data_all_out = sorted_data_all_out[,-3]
    }
    else{
      sorted_data_out = merge(ica[,1:numb], sorted_data[2:nrow(sorted_data), ], by=0)
      sorted_data_all_out = merge(ica[,1:numb], sorted_data_all[2:nrow(sorted_data_all),], by=0)
    }
    sorted_data_out = rbind.fill(sorted_data[1,], sorted_data_out)
    sorted_data_all_out = rbind.fill(sorted_data_all[1,], sorted_data_all_out)
    colnames(sorted_data_out) = gsub("Row.names", nm, colnames(sorted_data_out))
    colnames(sorted_data_all_out) = gsub("Row.names", nm, colnames(sorted_data_all_out))
    write.csv(sorted_data_out, paste(outdir, n, "_", m, "_lite.csv", sep=""), row.names=FALSE)
    write.csv(sorted_data_all_out, paste(outdir, n, "_", m, ".csv", sep=""), row.names=FALSE)
  }
}


# proteomics
todolist("~/documents/Lund_Melanoma/Results/proteomics/ICA/ICA_proteomics_IC_Clinical_Correlation_P_Value_all.tsv",
         "~/documents/Lund_Melanoma/Results/proteomics/ICA/ICA_proteomics_IC_mean_mixing_score.txt",
         "~/documents/Lund_Melanoma/Results/proteomics/ICA/significant_IC_clinical.csv")
GSEA("~/documents/Lund_Melanoma/Results/proteomics/ICA/MG_ICA_proteomics_IC_centroid.csv", 
     "~/documents/Lund_Melanoma/Results/proteomics/ICA/significant_IC_clinical.csv", 
     "~/documents/Lund_Melanoma/Results/proteomics/GSEA/")
HMP(centroid_file = "~/documents/Lund_Melanoma/Results/proteomics/ICA/MG_ICA_proteomics_IC_centroid.csv",
    clinical_file = "~/documents/Lund_Melanoma/Data/proteomics/proteomics_clinical.csv",
    data_file = "~/documents/Lund_Melanoma/Data/proteomics/ICA_proteomics.csv",
    td_list_file = "~/documents/Lund_Melanoma/Results/proteomics/ICA/significant_IC_clinical.csv",
    outdir = "~/documents/Lund_Melanoma/Results/proteomics/GSEA/", numb = 2, nm = "Accession")

# transcriptomics
todolist("~/documents/Lund_Melanoma/Results/transcriptomics/ICA/ICA_transcriptomics_IC_Clinical_Correlation_P_Value_all.tsv",
         "~/documents/Lund_Melanoma/Results/transcriptomics/ICA/ICA_transcriptomics_IC_mean_mixing_score.txt",
         "~/documents/Lund_Melanoma/Results/transcriptomics/ICA/significant_IC_clinical.csv")
GSEA("~/documents/Lund_Melanoma/Results/transcriptomics/ICA/MG_ICA_transcriptomics_IC_centroid.csv", 
     "~/documents/Lund_Melanoma/Results/transcriptomics/ICA/significant_IC_clinical.csv", 
     "~/documents/Lund_Melanoma/Results/transcriptomics/GSEA/")
HMP(centroid_file = "~/documents/Lund_Melanoma/Results/transcriptomics/ICA/MG_ICA_transcriptomics_IC_centroid.csv",
    clinical_file = "~/documents/Lund_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv",
    data_file = "~/documents/Lund_Melanoma/Data/transcriptomics/ICA_transcriptomics.csv",
    td_list_file = "~/documents/Lund_Melanoma/Results/transcriptomics/ICA/significant_IC_clinical.csv",
    outdir = "~/documents/Lund_Melanoma/Results/transcriptomics/GSEA/", numb = 1, nm="Gene.name")


# phospho
todolist("~/documents/Lund_Melanoma/Results/phospho/ICA/ICA_phospho_IC_Clinical_Correlation_P_Value_all.tsv",
         "~/documents/Lund_Melanoma/Results/phospho/ICA/ICA_phospho_IC_mean_mixing_score.txt",
         "~/documents/Lund_Melanoma/Results/phospho/ICA/significant_IC_clinical.csv")
GSEA("~/documents/Lund_Melanoma/Results/phospho/ICA/MG_ICA_phospho_IC_centroid.csv", 
     "~/documents/Lund_Melanoma/Results/phospho/ICA/significant_IC_clinical.csv", 
     "~/documents/Lund_Melanoma/Results/phospho/GSEA/")
HMP(centroid_file = "~/documents/Lund_Melanoma/Results/phospho/ICA/MG_ICA_phospho_IC_centroid.csv",
    clinical_file = "~/documents/Lund_Melanoma/Data/phospho/phospho_clinical.csv",
    data_file = "~/documents/Lund_Melanoma/Data/phospho/ICA_phospho.csv",
    td_list_file = "~/documents/Lund_Melanoma/Results/phospho/ICA/significant_IC_clinical.csv",
    outdir = "~/documents/Lund_Melanoma/Results/phospho/GSEA/", numb = 3, nm="Modified_sequence")





