# KM check PIK3CB phosphosites
library(dplyr)
library(readr)

phospho <- read_tsv("Data/phospho/MM_Phosho_DIA_lund_Cohort_total_peptides_NORMALIZED.txt")
survival <- read_csv("Data/clinical&histology.csv")


survival <- survival %>%
  select(samples=X1, dss.days, survival_5yr,	survival_3yr,	survival_1yr, survival_6mo)

phospho.filt <- phospho %>%
  filter(grepl("PIK3C", Gene)) %>%
  select(-c("PEP.StrippedSequence", "PG.ProteinAccessions", "PG.ProteinDescriptions", "PG.ProteinGroups", "PG.ProteinNames", "PG.Qvalue"))

survival <- survival %>%
  filter(samples %in% intersect(colnames(phospho.filt[3:124]), survival$samples))

phospho.filt <- phospho.filt[c("EG.ModifiedSequence", "Gene", intersect(colnames(phospho.filt[3:124]), survival$samples))]
phospho.filt["phosphosites"] = paste(phospho.filt$Gene, phospho.filt$EG.ModifiedSequence, sep='')
phospho.filt = phospho.filt[, c(118, 3:117)]
phospho.filt = t(phospho.filt)
colnames(phospho.filt) = phospho.filt[1, ]
phospho.filt = phospho.filt[-1, ]
phospho.filt = data.frame(phospho.filt)
phospho.filt$samples = unlist(rownames(phospho.filt))

survival_phos <- survival %>%
  inner_join(phospho.filt, by="samples")

write.csv(survival_phos, "Data/phospho/PIK3C_survival.csv", row.names=FALSE)




#====INSTALL PACKAGES======================================

# List of packages to install

.packages = c("devtools","ggplot2","ggbiplot","RColorBrewer","bitops",
              "dplyr","ggpubr","tidyverse",
              "survival","prodlim","ranger", "ggfortify","survminer",
              "pROC","nlme","reshape","dplyr","FactoMineR","verification","pastecs")

# Install packages if they are not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  install_github("vqv/ggbiplot")}

# load the packages if they are not loaded
lapply(.packages, require, character.only=TRUE)
#=================================================================

survival_phos <- read_csv("Data/phospho/PIK3C_survival.csv")
survival_phos <- survival_phos[, which(colMeans(!is.na(survival_phos)) > 0.3)]
survival_phos[, 3:6] <- abs(survival_phos[, 3:6]-1) 
roc.5yr <- survival_phos[, c(1:3, 7:14)]
roc.5yr <- roc.5yr[!is.na(roc.5yr[,3]),]
roc.3yr <- survival_phos[, c(1,2,4, 7:14)]
roc.3yr <- roc.3yr[!is.na(roc.3yr[,3]),]
roc.1yr <- survival_phos[, c(1,2,5, 7:14)]
roc.1yr <- roc.1yr[!is.na(roc.1yr[,3]),]
roc.6mo <- survival_phos[, c(1,2,6, 7:14)]
roc.6mo <- roc.6mo[!is.na(roc.6mo[,3]),]


km.roc <- function(DATA, cutoff.cat.roc){
  #=== ROC analysis
  
  colnames1 <- as.character(colnames(DATA))
  ROC.matrix <- matrix(nrow=ncol(DATA)-3,ncol = 6)
  for(i in 4:ncol(DATA)) {
    DATA.X = na.omit(DATA[, c(1:3, i)])
    
    png(filename = paste("Results/phospho/PIK3C_ROC/",cutoff.cat.roc,"/roc_",colnames1[i],".png",sep=""), width = 500, height = 480, units = "px", pointsize = 16,
        bg = "white")
    
    
    roc_i <- roc(unlist(DATA.X[,3]) ~ unlist(DATA.X[,4]),data=DATA.X,percent=T,
                 smooth=F,
                 legacy.axes=F,
                 ci=F, 
                 boot.n=2000, 
                 ci.alpha=0.7, 
                 stratified=T,
                 plot=TRUE, 
                 auc.polygon=T, 
                 max.auc.polygon=T, 
                 grid=TRUE,
                 print.auc=TRUE, 
                 print.thres='best', 
                 print.thres.best.method='y',   # "closest.topleft", "youden"
                 print.thres.adj=c(-0.7, 1.5),#,c(-0.05, 2.0),c(-0.05, 5.0)),
                 print.thres.pattern="Cut-off: %.2f    \n\n\nSpec: %.1f \nSens: %.1f",
                 print.thres.pattern.cex = 1.1,
                 print.auc.cex = 1.1,
                 print.auc.y=7,
                 print.auc.x=83,
                 cex.axis=1.5,
                 cex.lab=1.5,
                 print.thres.pch=16,
                 print.thres.cex=2.0,
                 cex.main=1.5,
                 na.rm=T)
    
    print(roc_i)
    dev.off()
    
    ROC.matrix[i-3,1]<- roc_i$auc
    
    ROC_table <- cbind(roc_i$thresholds, roc_i$specificities, roc_i$sensitivities)
    cut_off <- ROC_table[which.max(ROC_table[, 2] + ROC_table[, 3]), ] 
    
    roc.A <- verification::roc.area(roc_i$original.response,roc_i$original.predictor)
    
    ROC.matrix[i-3,2] <- roc.A$n.total  ##N number of samples
    ROC.matrix[i-3,3] <- roc.A$p.value  ##N p.value
    ROC.matrix[i-3,4] <- cut_off[1]     ##cut_off = best threshold based on "youden" method
    ROC.matrix[i-3,5] <- cut_off[2]     ##Sp
    ROC.matrix[i-3,6] <- cut_off[3]     ##Se
    
    
    roc_i.newlist <- list(AUC = as.matrix(roc_i$auc),
                          sens.= as.matrix(roc_i$sensitivities),
                          Spec.= as.matrix(roc_i$specificities),
                          Thesholds = as.matrix(roc_i$thresholds),
                          predictor = as.matrix(roc_i$predictor),
                          response = as.data.frame(roc_i$response))
    
    write.infile(roc_i.newlist, paste("Results/phospho/PIK3C_ROC/",cutoff.cat.roc,"/roc_",colnames1[i],".csv",sep=''), sep = "\t")
    
  }
  row.names(ROC.matrix) <- colnames(DATA)[4:ncol(DATA)]
  colnames(ROC.matrix) <- c("AUC(%)", "N","p.value","cut_off", "Sp", "Se")
  
  ROC.matrix <- as.data.frame(ROC.matrix)
  
  write.csv(ROC.matrix, paste("Results/phospho/PIK3C_ROC/",cutoff.cat.roc,"/ROC_",cutoff.cat.roc,".csv",sep=""))
  
  
  ### Kaplan meier
  
  km.matrix <- matrix(nrow=ncol(DATA)-3, ncol = 1, dimnames = list(colnames(DATA)[4:ncol(DATA)],c("pv.km")))
  
  for (j in 4:ncol(DATA)){
    
    data <- DATA[, c(1:3, j)]
    
    data1 <- na.omit(data)
    colnames(data1) <- c("samples", "time", "status", "Prot")
    
    cat.roc <- vector(mode = "numeric")
    
    for (i in 1:nrow(data1)) {
      if(data1[i, 4] < ROC.matrix$cut_off[j-3]) {cat.roc[i] <- "Low"} else {cat.roc[i] <- "High"}   # Low = 0, High = 1,  
    }
    
    data1$cat.roc <- cat.roc
    data1$cat.roc <- as.factor(data1$cat.roc)
    
    #KM
    
    km_fit <- survfit(Surv(time, status) ~ cat.roc, data=data1)
    
    pv <- survminer::surv_pvalue(km_fit, method = "log-rank", test.for.trend = F, combine = F, data = data1)      #https://urldefense.proofpoint.com/v2/url?u=https-3A__rpkgs.datanovia.com_survminer_reference_surv-5Fpvalue.html&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=vdmEsUtCQDmg1Lp_kZJV_VVeYj4jchLUD9dty404uE4&m=-_wDTGk9iAFH5ZoI6wMcQeTbgV0hLzaVJk_fVKHZClA&s=1jpjzrsgB6HIUNKasLY89uCUPGL84EGV-hf7NoEoyIQ&e= 
    
    km.matrix[j-3,1] <- pv$pval
    
    
    png(filename = paste("Results/phospho/PIK3C_KM/",cutoff.cat.roc,"/km_",colnames(data)[4],".png",sep =""), width = 800, height = 500, units = "px", pointsize = 16,
        bg = "white")
    
    km.plot <- autoplot(km_fit, conf.int = F, censor = T, surv.size = 1)+theme_bw()+
      labs(title = paste("km_",colnames(data)[4], "  p=",signif(pv$pval, 4),
                         " |  auc=",signif(ROC.matrix$`AUC(%)`[j], 3),"%",
                         "  cut-off=",signif(ROC.matrix$cut_off[j], 2),sep =""),
           y = paste("Survival (",cutoff.cat.roc,")",sep=""), 
           x = "Time (days)")+
      theme (text = element_text(size=13), legend.title = element_blank())
    #km.plot
    print(km.plot)
    
    dev.off()
  }
  km.matrix <- as.data.frame(km.matrix)
  km.matrix$sig <- cut(km.matrix[,1], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  
  write.csv(km.matrix, paste("Results/phospho/PIK3C_KM/",cutoff.cat.roc,"/Km_",cutoff.cat.roc,".csv",sep=""))
  return(km.matrix)
}

km.results <- km.roc(DATA=roc.5yr, cutoff.cat.roc = "5yr")
km.results <- km.roc(DATA=roc.3yr, cutoff.cat.roc = "3yr")
km.results <- km.roc(DATA=roc.1yr, cutoff.cat.roc = "1yr")




