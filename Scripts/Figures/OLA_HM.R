library(dplyr)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Heatmap for OLA proteomics
OLA = read.csv("~/Documents/Segundo_Melanoma/Results/OLA_summary.csv")
prot.clinical = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names = 1)
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")

prot.clinical = prot.clinical[order(prot.clinical['dss.days']), ]
prot.clinical.1 = prot.clinical[, c(2, 13, 23)]

OLA.prot = OLA[OLA['Group'] == "proteomics", -c(2,4,6)]
OLA.prot$Enriched_in = paste(OLA.prot$Feature, OLA.prot$Enriched_in, sep='_')
OLA.prot = OLA.prot[, -4]
prot.1 = prot[, -112]
m = match(rownames(prot.clinical.1), colnames(prot.1),)
m[111] = 111
m[112] = 112
prot.1 = prot.1[, m]

OLA.prot.table <- merge(OLA.prot,prot.1,by=c("Gene.name","Accession"))
OLA.prot.table = OLA.prot.table[order(OLA.prot.table$Enriched_in), ]
OLA.prot.table$FDR = round(OLA.prot.table$FDR, 5) 

OLA.prot.table.out = rbind.fill(OLA.prot.table, as.data.frame(t(prot.clinical.1)))
OLA.prot.table.out = OLA.prot.table.out[c(c(95:97), c(1:94)), ]
row.names(OLA.prot.table.out)[1:3] = c('stage', 'collection_to_death/censor', 'NRAS')
write.csv(OLA.prot.table.out, "~/documents/Segundo_Melanoma/Results/OLA_HM_proteomics.csv", row.names = TRUE)

col_fun1 = colorRamp2(c(0, 0.05), c("white", "green"))
col_fun2 = colorRamp2(c(0, 10000), c("white","blue"))
col_fun3 = colorRamp2(c(1,4), c("white", "orange"))
col_fun4 = colorRamp2(c(0, 1), c("white", "purple"))
gn = rowAnnotation(FDR = OLA.prot.table$FDR,
                   gene.name = anno_text(OLA.prot.table$Gene.name,
                                         location = 0.5, just = "center"),
                   accession = anno_text(OLA.prot.table$Accession,
                                         location = 0.5, just = "center"),
                   enriched.in = anno_text(OLA.prot.table$Enriched_in,
                                           location = 0.5, just = "center"),
                   col = list(FDR=col_fun1))

OLA.prot.table.1 = OLA.prot.table %>% select(matches("MM"))
breaksList = seq(min(prot.1[,1:110]), max(prot.1[,1:110]), by=1)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
anno = HeatmapAnnotation(days = as.numeric(prot.clinical.1[,2]), 
                       stage = as.numeric(prot.clinical.1[,1]), 
                       NRAS = as.numeric(prot.clinical.1[,3]),
                       col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

pdf("~/documents/Segundo_Melanoma/Results/OLA_HM_proteomics.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(OLA.prot.table.1), col = col, column_title = paste("proteomics outliers"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
             cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

pdf("~/documents/Segundo_Melanoma/Results/both_OLA_HM_proteomics.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(OLA.prot.table.1), col = col, column_title = paste("proteomics outliers"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

pdf("~/documents/Segundo_Melanoma/Results/row_OLA_HM_proteomics.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(OLA.prot.table.1), col = col, column_title = paste("proteomics outliers"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
             cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

pdf("~/documents/Segundo_Melanoma/Results/col_OLA_HM_proteomics.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(OLA.prot.table.1), col = col, column_title = paste("proteomics outliers"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
             cluster_rows = FALSE, cluster_columns = TRUE, show_row_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()


# Heatmap for OLA transcriptomics
trans.clinical = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names = 1)
OLA = read.csv("~/Documents/Segundo_Melanoma/Results/OLA_summary.csv")
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")

trans.clinical = trans.clinical[order(trans.clinical['dss.days']), ]
trans.clinical.1 = trans.clinical[, c(2, 13, 23)]

OLA.trans = OLA[OLA['Group'] == "transcriptomics", -c(2,3,4,6)]
OLA.trans$Enriched_in = paste(OLA.trans$Feature, OLA.trans$Enriched_in, sep='_')
OLA.trans = OLA.trans[, -3]
trans.1 = trans[, -2]
m = match(rownames(trans.clinical.1), colnames(trans.1),)
m[135] = 1
trans.1 = trans.1[, m]

OLA.trans.table <- merge(OLA.trans,trans.1,by=c("Gene.name"))
OLA.trans.table = OLA.trans.table[order(OLA.trans.table$Enriched_in), ]
OLA.trans.table$FDR = round(OLA.trans.table$FDR, 5) 

OLA.trans.table.out = rbind.fill(OLA.trans.table, as.data.frame(t(trans.clinical.1)))
OLA.trans.table.out = OLA.trans.table.out[c(c(21:23), c(1:20)), ]
row.names(OLA.trans.table.out)[1:3] = c('stage', 'collection_to_death/censor', 'NRAS')
write.csv(OLA.trans.table.out, "~/documents/Segundo_Melanoma/Results/OLA_HM_transcriptomics.csv", row.names = TRUE)

col_fun1 = colorRamp2(c(0, 0.05), c("white", "green"))
col_fun2 = colorRamp2(c(0, 10000), c("white","blue"))
col_fun3 = colorRamp2(c(1,4), c("white", "orange"))
col_fun4 = colorRamp2(c(0, 1), c("white", "purple"))
gn = rowAnnotation(FDR = OLA.trans.table$FDR,
                   gene.name = anno_text(OLA.trans.table$Gene.name,
                                         location = 0.5, just = "center"),
                   enriched.in = anno_text(OLA.trans.table$Enriched_in,
                                           location = 0.5, just = "center"),
                   col = list(FDR=col_fun1))

OLA.trans.table.1 = OLA.trans.table %>% select(matches("MM"))
breaksList = seq(min(trans.1[,1:134]), max(trans.1[,1:134]), by=1)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
anno = HeatmapAnnotation(days = as.numeric(trans.clinical.1[,2]), 
                         stage = as.numeric(trans.clinical.1[,1]), 
                         NRAS = as.numeric(trans.clinical.1[,3]),
                         col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

pdf("~/documents/Segundo_Melanoma/Results/OLA_HM_transcriptomics.pdf", height = 6.5, width = 30)
hp = Heatmap(as.matrix(OLA.trans.table.1), col = col, column_title = paste("transcriptomics outliers"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
             cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

# Heatmap for OLA phospho
phospho.clinical = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names = 1)
OLA = read.csv("~/Documents/Segundo_Melanoma/Results/OLA_summary.csv")
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")

phospho.clinical = phospho.clinical[order(phospho.clinical['dss.days']), ]
phospho.clinical.1 = phospho.clinical[, c(2, 13, 23)]

OLA.phospho = OLA[OLA['Group'] == "phospho", -c(2,4)]
OLA.phospho$Enriched_in = paste(OLA.phospho$Feature, OLA.phospho$Enriched_in, sep='_')
OLA.phospho = OLA.phospho[, -5]
phospho.1 = phospho[, -3]
m = match(rownames(phospho.clinical.1), colnames(phospho.1),)
m[119] = 1
m[120] = 2
m[121] = 3
phospho.1 = phospho.1[, m]

OLA.phospho.table <- merge(OLA.phospho,phospho.1,by=c("Modified_sequence", "Gene.name", "Accession"))
OLA.phospho.table = OLA.phospho.table[order(OLA.phospho.table$Enriched_in), ]
OLA.phospho.table$FDR = round(OLA.phospho.table$FDR, 5) 

OLA.phospho.table.out = rbind.fill(OLA.phospho.table, as.data.frame(t(phospho.clinical.1)))
OLA.phospho.table.out = OLA.phospho.table.out[c(c(5:7), c(1:4)), ]
row.names(OLA.phospho.table.out)[1:3] = c('stage', 'collection_to_death/censor', 'NRAS')
write.csv(OLA.phospho.table.out, "~/documents/Segundo_Melanoma/Results/OLA_HM_phospho.csv", row.names = TRUE)

col_fun1 = colorRamp2(c(0, 0.05), c("white", "green"))
col_fun2 = colorRamp2(c(0, 10000), c("white","blue"))
col_fun3 = colorRamp2(c(1,4), c("white", "orange"))
col_fun4 = colorRamp2(c(0, 1), c("white", "purple"))
gn = rowAnnotation(FDR = OLA.phospho.table$FDR,
                   gene.name = anno_text(OLA.phospho.table$Gene.name,
                                         location = 0.5, just = "center"),
                   accession = anno_text(OLA.phospho.table$Accession,
                                         location = 0.5, just = "center"),
                   enriched.in = anno_text(OLA.phospho.table$Enriched_in,
                                           location = 0.5, just = "center"),
                   sequence = anno_text(OLA.phospho.table$Modified_sequence,
                                        location = 0.5, just = "center"),
                   col = list(FDR=col_fun1))

OLA.phospho.table.1 = OLA.phospho.table %>% select(matches("MM"))
breaksList = seq(min(phospho.1[,1:118]), max(phospho.1[,1:118]), by=1)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
anno = HeatmapAnnotation(days = as.numeric(phospho.clinical.1[,2]), 
                         stage = as.numeric(phospho.clinical.1[,1]), 
                         NRAS = as.numeric(phospho.clinical.1[,3]),
                         col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

pdf("~/documents/Segundo_Melanoma/Results/OLA_HM_phospho.pdf", height = 4, width = 30)
hp = Heatmap(as.matrix(OLA.phospho.table.1), col = col, column_title = paste("phosphoeomics outliers"), top_annotation = anno,  right_annotation=gn, show_column_names = FALSE,
             cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()
