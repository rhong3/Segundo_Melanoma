
##https://urldefense.proofpoint.com/v2/url?u=https-3A__rviews.rstudio.com_2017_09_25_survival-2Danalysis-2Dwith-2Dr_&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=vdmEsUtCQDmg1Lp_kZJV_VVeYj4jchLUD9dty404uE4&m=-_wDTGk9iAFH5ZoI6wMcQeTbgV0hLzaVJk_fVKHZClA&s=LUaiy7QogkhUGGQHjylGXPsj6QT09q8ePdlfhrU8vAs&e= 

setwd("~/R/Rstudio/Segundo/Kaplan-Meier analysis ROC curve/test2") 


Prot.data <- read.delim("20200716_ROC_KM_phospho_top_list_candidates.txt", header = TRUE, row.names = 1)
head(Prot.data, 5)



Prot.data1 <- Prot.data[, 2:ncol(Prot.data)] 



#new update 2020-05-20
covariates1 <- read.delim("20200520ClinicalData_144samplesLundMM_REVISEDandUPDATEDtill2020_no_formula.txt", header = TRUE,row.names = 1)
covariates1$sample <- row.names(covariates1)

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

roc.Gr <- as.data.frame(covariates1[, c("DSS.DAYS.v1.days.from.first.metastasis.to.death.or.censoring","dss.events","sample")])   # dataframe to use in the ROC analysis

roc.Gr <- na.omit(roc.Gr)  # remove missing values

####___Parameters of the function "km.roc"

#               m = dataframe with expresion data (e.g proteomcis data; columns are samples and rows are proteins)
#          groups = dataframe with survival covariates (e.g. 'roc.Gr', columns = c("DSS.DAYS.","dss.events","sample")). 
#                   Must include a column called "sample" containing the sample codes.
#            gene = vector with the list of gene code per proteins
#    time.nameCol = name of the column that indicates the survival time (e.g. "OS..Days")
#  status.nameCol = name of the column that indicates the censored status (dead (1) or not dead (0))(e.g. "Status.Dead.alive"). 
#  cutoff.cat.roc = Time in years or month in which the samples shoud be separated in two groups 
#                   (example: if with want to separate the samples into <3 years survival and >3 years survival, the cutoff.cat.roc = 3)
#     time.format = according to the cutoff.cat.roc, 'time.format' indicates if the time is given in "years" or "months".
#                   (example: time.format = "years")


km.roc <- function(m,groups,gene,time.nameCol,status.nameCol,cutoff.cat.roc,time.format){
  
  dir.create("ROC results1")
  dir.create("Km_graph1")
  
  roc.Gr1 <- groups[groups[,status.nameCol] %in% 1,]   ### only deads
  
  
  # create categorical survival columns (e.g. >3y or <3y) 
  
  if(time.format == "years") d=365.35 else d=30.44
  
  roc.Gr1$time.f <- roc.Gr1[,time.nameCol]/d
  
  cat.surv <- vector(mode = "numeric")
  
  for (i in 1:nrow(roc.Gr1)) {
    if(roc.Gr1$time.f[i] < cutoff.cat.roc) cat.surv[i] <- 1 else cat.surv[i] <- 0
  }
  
  roc.Gr1$cat.surv <- cat.surv
  #roc.Gr1$cat.surv <- as.factor(roc.Gr1$cat.surv)
  
  m1 <- as.data.frame(t(m))
  m1$sample <- row.names(m1)
  
  DATA <- plyr::join_all(list(roc.Gr1[,c("sample", "cat.surv")],m1), by="sample")   # Merging data.frames (example: "m1" AND "roc.Gr1[,c("sample", "cat.surv")]")
  row.names(DATA) <- DATA$sample
  DATA <- DATA[,-1]
  
  #=== ROC analysis
  
    colnames1 <- as.character(colnames(DATA))
    ROC.matrix <- matrix(nrow=ncol(DATA),ncol = 6)
    p= 0
    for(i in 2:ncol(DATA)) {
      
      
      png(filename = paste("ROC results1/","roc_",colnames1[i],"_",gene[i-1],".png",sep=""),width = 500, height = 480, units = "px", pointsize = 16,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
      
      
      roc_i <- roc(DATA[,1] ~ DATA[,i],data=DATA,percent=T,
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
      
      ROC.matrix[i,1]<- roc_i$auc
      
      ROC_table <- cbind(roc_i$thresholds, roc_i$specificities, roc_i$sensitivities)
      cut_off <- ROC_table[which.max(ROC_table[, 2] + ROC_table[, 3]), ] 
      
      roc.A <- verification::roc.area(roc_i$original.response,roc_i$original.predictor)
      
      ROC.matrix[i,2] <- roc.A$n.total  ##N number of samples
      ROC.matrix[i,3] <- roc.A$p.value  ##N p.value
      ROC.matrix[i,4] <- cut_off[1]     ##cut_off = best threshold based on "youden" method
      ROC.matrix[i,5] <- cut_off[2]     ##Sp
      ROC.matrix[i,6] <- cut_off[3]     ##Se
      
      
      roc_i.newlist <- list(AUC = as.matrix(roc_i$auc),
                            sens.= as.matrix(roc_i$sensitivities),
                            Spec.= as.matrix(roc_i$specificities),
                            Thesholds = as.matrix(roc_i$thresholds),
                            predictor = as.matrix(roc_i$predictor),
                            response = as.data.frame(roc_i$response))
      
      #write.infile(roc_i.newlist, paste("ROC results1/","roc_",colnames1[i],"_",gene[i-1],".csv",sep=''), sep = "\t")
      
      p= p+1
    }
    row.names(ROC.matrix) <- colnames(DATA)
    colnames(ROC.matrix) <- c("AUC(%)", "N","p.value","cut_off", "Sp", "Se")
    
    ROC.matrix <- as.data.frame(ROC.matrix[-1,])
    ROC.matrix$gene <- gene
    
    write.csv(ROC.matrix, paste("ROC results1/","ROC_",cutoff.cat.roc," ",time.format,".csv",sep=""))
    
  
  ### Kaplan meier
    
  km.matrix <- matrix(nrow=nrow(m),ncol = 1, dimnames = list(row.names(m),c("pv.km")))
  
  for (j in 1:nrow(m)){
  
          prot.data <- as.data.frame(t(m[j,]))
          prot.data$sample <- row.names(prot.data)
          
          
          data <- plyr::join_all(list(prot.data,groups),by="sample")
          row.names(data) <- data$sample
          data <- data[,-2]
                   
          data1 <- na.omit(data)
          colnames(data1) <- c("Prot","time", "status")
          
          cat.roc <- vector(mode = "numeric")
              
          for (i in 1:nrow(data1)) {
              if(data1[i, 1] < ROC.matrix$cut_off[j]) {cat.roc[i] <- "Low"} else {cat.roc[i] <- "High"}   # Low = 0, High = 1,  
            }
              
          data1$cat.roc <- cat.roc
          data1$cat.roc <- as.factor(data1$cat.roc)
              
          data1$sample <- row.names(data1)
              
          #KM
              
          km_fit <- survfit(Surv(time, status) ~ cat.roc, data=data1)
              
          pv <- survminer::surv_pvalue(km_fit, method = "log-rank", test.for.trend = F, combine = F, data = data1)      #https://urldefense.proofpoint.com/v2/url?u=https-3A__rpkgs.datanovia.com_survminer_reference_surv-5Fpvalue.html&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=vdmEsUtCQDmg1Lp_kZJV_VVeYj4jchLUD9dty404uE4&m=-_wDTGk9iAFH5ZoI6wMcQeTbgV0hLzaVJk_fVKHZClA&s=1jpjzrsgB6HIUNKasLY89uCUPGL84EGV-hf7NoEoyIQ&e= 
                                
          km.matrix[j,1] <- pv$pval
          
              
          png(filename = paste("Km_graph1/","km_",colnames(data)[1],"_",gene[j],".png",sep =""), width = 500, height = 400, units = "px", pointsize = 16,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
              
          # ggsurvplot(km_fit)
          
          km.plot <- autoplot(km_fit, conf.int = F, censor = T, surv.size = 1)+theme_bw()+
                     labs(title = paste("km_",colnames(data)[1],"_",gene[j], "  p=",signif(pv$pval, 4),
                                        " |  auc=",signif(ROC.matrix$`AUC(%)`[j], 3),"%",
                                        "  cut-off=",signif(ROC.matrix$cut_off[j], 2),sep =""),
                          y = paste("Survival (",cutoff.cat.roc," ",time.format,")",sep=""), 
                          x = "Time (days)")+
                      theme (text = element_text(size=13), legend.title = element_blank())
          #km.plot
          print(km.plot)
               
          dev.off()
   }
  km.matrix <- as.data.frame(km.matrix)
  km.matrix$gene <- gene
  km.matrix$sig <- cut(km.matrix[,1], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  
  write.csv(km.matrix, paste("Km_graph1/","Km_",cutoff.cat.roc," ",time.format,".csv",sep=""))
  return(km.matrix)
}

km.results <- km.roc(m=Prot.data1, 
                     groups=roc.Gr,
                     gene=Prot.data$Gene.name, 
                     time.nameCol="DSS.DAYS.v1.days.from.first.metastasis.to.death.or.censoring",
                     status.nameCol="dss.events", 
                     cutoff.cat.roc = 5,
                     time.format = "years")   # "years" or "months".
