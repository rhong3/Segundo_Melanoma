# ICA and proteomics subtypes
library(readxl)
library(fastDummies)
library(dplyr)

clinical = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv")
subtype = read_xlsx("~/documents/Segundo_Melanoma/Data/Sample name and subtypes.xlsx")
colnames(clinical)[1] = 'Sample name'

subtype = dummy_cols(subtype,  select_columns = 'SUBTYPE')

joined = left_join(clinical, subtype, by='Sample name')[, c(1, 52:56)]
joined[is.na(joined)] = 0
row.names(joined) = joined$`Sample name`
colnames(joined) = gsub("SUBTYPE_", "", colnames(joined))
joined = joined[,-1]
write.csv(joined, '~/documents/Segundo_Melanoma/Data/prot_subtypes.csv')
