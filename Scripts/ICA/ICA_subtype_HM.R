# ICA subtype Heatmap
library(BiocManager)
library('org.Hs.eg.db')
library("ReactomePA")
library("fgsea")
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(circlize)
library("RColorBrewer")
library(colorspace)


centroid_file = "~/documents/Segundo_Melanoma/Results/proteomics/ICA/MG_ICA_proteomics_IC_centroid.csv"
clinical_file = "~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv"
data_file = "~/documents/Segundo_Melanoma/Data/proteomics/ICA_proteomics.csv"
td_list_file = "~/documents/Segundo_Melanoma/Results/proteomics/ICA/0.00001/significant_IC_clinical.csv"
outdir = "~/documents/Segundo_Melanoma/Results/proteomics/GSEA/strict/" 
numb = 2 
nm = "Accession"

ica = read.csv(centroid_file, row.names=1)
clinical = read.csv(clinical_file, row.names=1)
clinical$prim.breslow.class = gsub("NM", "NA", clinical$prim.breslow.class)
colnames(clinical) = gsub("lymphnode..", "lymphnode.percentage", colnames(clinical))
colnames(clinical) = gsub("tumor..", "tumor.percentage", colnames(clinical))
colnames(clinical) = gsub("tumor.percentageell.size.average", "tumor.cell.size.average", colnames(clinical))
ori_data = read.csv(data_file, row.names=1)
ori_data[ori_data>5] <- 5
ori_data[ori_data<(-5)] <- -5
td_list <- read.csv(td_list_file)
td_list = td_list[td_list$Feature %in% c('EC.Mit'),]
for (a in 1:nrow(td_list)){
  m = toString(droplevels(td_list[a, "Feature"]))
  n = toString(droplevels(td_list[a, "IC"]))
  ica = ica[order(-ica[n]), ]
  clinical = clinical[order(clinical[m]), ]
  sorted_data_all = ori_data[match(rownames(ica), rownames(ori_data)), match(rownames(clinical), colnames(ori_data))]
  sorted_data = data.matrix(sorted_data_all[-c(11:(nrow(sorted_data_all))),])
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
  col = colorRampPalette(colorspace::diverge_hsv(10))(11)[breaksList+6]
  col_fun = c("0" = "white", "1" = "#DECBE4")
  anno = HeatmapAnnotation('EC.Mit' = as.factor(sorted_data_out[1, ]), col=list('EC.Mit'=col_fun), annotation_legend_param = list(direction = "horizontal"))
  pdf(paste(outdir, n, "_", m, "_HM.pdf", sep=""), height = 3, width = 20)
  hp = Heatmap(as.matrix(sorted_data_out[2:nrow(sorted_data_out), ]), col = col, column_title = paste(n), 
               column_title_gp = gpar(fontsize = 18, fontface = "bold"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
               cluster_rows = FALSE, cluster_columns = FALSE, name = "Value", heatmap_legend_param = list(direction = "horizontal"))
  draw(hp, heatmap_legend_side = "bottom", 
       annotation_legend_side = "bottom", merge_legend = TRUE,)
  dev.off()
}

