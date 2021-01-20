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



