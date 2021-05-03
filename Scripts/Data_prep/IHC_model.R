# IHC survival models
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(survival)

# IHC data filter
DDX11 <- read_excel("~/Documents/Segundo_Melanoma/Data/DDX11_IHC.xlsx")
DDX11 = DDX11[, c(1:6)]
DDX11.mean = DDX11 %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)


NBP1 <- read_excel("~/Documents/Segundo_Melanoma/Data/NBP1_IHC.xls")
NBP1 = NBP1[, c(1:7)]
NBP1.mean = NBP1 %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

ADAM10 <- read_excel("~/Documents/Segundo_Melanoma/Data/ADAM10_IHC.xls")
ADAM10 = ADAM10[, c(1:7)]
ADAM10.mean = ADAM10 %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

PIK3CB <- read_excel("~/Documents/Segundo_Melanoma/Data/PIK3CB_IHC.xls")
PIK3CB = PIK3CB[, c(1:7)]
PIK3CB.mean = PIK3CB %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

PAEP <- read_excel("~/Documents/Segundo_Melanoma/Data/PAEP_IHC.xls")
PAEP = PAEP[, c(1:7)]
PAEP.mean = PAEP %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

FGA <- read_excel("~/Documents/Segundo_Melanoma/Data/FGA_IHC.xls")
FGA = FGA[, c(1:7)]
FGA.mean = FGA %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

CDK4 <- read_excel("~/Documents/Segundo_Melanoma/Data/CDK4_IHC.xls")
CDK4 = CDK4[, c(1:7)]
CDK4.mean = CDK4 %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

HMOX1 <- read_excel("~/Documents/Segundo_Melanoma/Data/HMOX1_IHC.xlsx")
HMOX1 = HMOX1[, c(1:3, 5:8)]
HMOX1.mean = HMOX1 %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

CTNND1 <- read_excel("~/Documents/Segundo_Melanoma/Data/CTNND1_IHC.xls")
CTNND1 = CTNND1[, c(1:7)]
CTNND1.mean = CTNND1 %>%
  fill(Chip_p) %>%
  na.omit() %>%
  group_by(Chip_p) %>%
  summarise_all(mean)

IHC = ADAM10.mean %>%
  left_join(PIK3CB.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(PAEP.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(FGA.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(CDK4.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(CTNND1.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(HMOX1.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(NBP1.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  left_join(DDX11.mean, by=c('Chip_p', 'DFS', 'PFS', 'OS', 'Live')) %>%
  select(c(1:3, 8:22, 4:7))

write.csv(IHC, "~/Documents/Segundo_Melanoma/Data/IHC_agg.csv", row.names = FALSE)
  

# Proteomics data filter
gene = c('PIK3CB', 'PAEP', 'FGA', 'CDK4', 'ADAM10', 'HMOX1', 'CTNND1', 'SCAI', 'DDX11')
proteomics <- read.csv("~/Documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
proteomics = proteomics[which(proteomics$Gene.name %in% gene), ]
proteomics = proteomics[,c(113, 1:110)]
row.names(proteomics) = proteomics$Gene.name
proteomics = t(proteomics)
proteomics = proteomics[-1,]

clinical = read.csv("~/Documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv")
row.names(clinical) = clinical[,1]
clinical = clinical[, c("dss.days", "dss.events", "survival_5yr",	"survival_3yr",	"survival_1yr", "survival_6mo")]

prot.clinical = merge(data.frame(proteomics), data.frame(clinical), by=0)
prot.clinical = na.omit(prot.clinical)
write.csv(prot.clinical, "~/Documents/Segundo_Melanoma/Data/IHC_prot_clinical.csv", row.names = FALSE)

# phospho data filter
gene = c('PIK3CB', 'PAEP', 'FGA', 'CDK4', 'ADAM10', 'HMOX1', 'CTNND1', 'SCAI', 'DDX11')
phospho = read.csv("~/Documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho = phospho[which(phospho$Gene.name %in% gene), ]
phospho = phospho[, c(4:122)]
phospho = phospho %>%
  group_by(Gene.name) %>%
  summarise_all(mean)
row.names(phospho) = phospho$Gene.name
phospho = t(phospho)
phospho = phospho[-1,]

clinical = read.csv("~/Documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv")
row.names(clinical) = clinical[,1]
clinical = clinical[, c("dss.days", "dss.events", "survival_5yr",	"survival_3yr",	"survival_1yr", "survival_6mo")]

phos.clinical = merge(data.frame(phospho), data.frame(clinical), by=0)
phos.clinical = na.omit(phos.clinical)
write.csv(phos.clinical, "~/Documents/Segundo_Melanoma/Data/IHC_phos_clinical.csv", row.names = FALSE)


# LM 
# prep
prot.clinical = read.csv("~/Documents/Segundo_Melanoma/Data/IHC_prot_clinical.csv")
phos.clinical = read.csv("~/Documents/Segundo_Melanoma/Data/IHC_phos_clinical.csv")
IHC = read.csv("~/Documents/Segundo_Melanoma/Data/IHC_agg.csv")
colnames(IHC) = gsub("PIK3cB", "PIK3CB", colnames(IHC))
colnames(IHC) = gsub("HMOX", "HMOX1", colnames(IHC))
# gene_means = colMeans(prot.clinical[,2:9])
# gene_std = apply(prot.clinical[,2:9], 2, sd)
# IHC.scale = IHC
# center_factors = c(0.0302873, 0.0302873, 0.0896912, 0.0896912, -0.2300511, -0.2300511, -0.1936995, -0.1936995, 0.1471385, 0.1471385, -0.2320132, -0.2320132, 0.1701781, 0.1701781)
# scale_factors = c(0.3197894, 0.3197894, 0.4426608, 0.4426608, 1.4249570, 1.4249570, 1.3063025, 1.3063025, 0.5791502, 0.5791502, 0.5169386, 0.5169386, 0.6680108, 0.6680108)
# IHC.scale[, 2:17] = scale(IHC.scale[, 2:17], center=scale_factors, scale=FALSE)

# plot
IHC.plot = data.frame()
for (i in 2:18){
  IHC.plot.temp = data.frame(Values=IHC[,i], Measure=colnames(IHC)[i], Days=IHC$OS, Live=IHC$Live)
  IHC.plot = rbind.data.frame(IHC.plot, IHC.plot.temp)
}
IHC.plot$Live = as.factor(IHC.plot$Live)
ggplot(IHC.plot, aes(Days, Values, color=Live)) +
  geom_point() +
  facet_wrap(~ Measure, ncol=6) + ggtitle("IHC vs. Survival") + theme_bw()

# IHC.scale.plot = data.frame()
# for (i in 2:15){
#   IHC.scale.plot.temp = data.frame(Values=IHC.scale[,i], Measure=colnames(IHC.scale)[i], Days=IHC.scale$OS, Live=IHC.scale$Live)
#   IHC.scale.plot = rbind.data.frame(IHC.scale.plot, IHC.scale.plot.temp)
# }
# IHC.scale.plot$Live = as.factor(IHC.scale.plot$Live)
# ggplot(IHC.scale.plot, aes(Days, Values, color=Live)) +
#   geom_point() +
#   facet_wrap(~ Measure, ncol=7) + ggtitle("IHC-scaled vs. Survival") + theme_bw()

prot.clinical.plot=data.frame()
for (i in 2:10){
  prot.clinical.temp = data.frame(Values=prot.clinical[,i], Measure=colnames(prot.clinical)[i], Days=prot.clinical$dss.days, Live=abs(prot.clinical$dss.events-1))
  prot.clinical.plot = rbind.data.frame(prot.clinical.plot, prot.clinical.temp)
}
prot.clinical.plot$Live = as.factor(prot.clinical.plot$Live)
ggplot(prot.clinical.plot, aes(Days, Values, color=Live)) +
  geom_point() +
  facet_wrap(~ Measure, ncol=3) + ggtitle("Proteome vs. Survival") + theme_bw()

phos.clinical.plot=data.frame()
for (i in 2:5){
  phos.clinical.temp = data.frame(Values=phos.clinical[,i], Measure=colnames(phos.clinical)[i], Days=phos.clinical$dss.days, Live=abs(phos.clinical$dss.events-1))
  phos.clinical.plot = rbind.data.frame(phos.clinical.plot, phos.clinical.temp)
}
phos.clinical.plot$Live = as.factor(phos.clinical.plot$Live)
ggplot(phos.clinical.plot, aes(Days, Values, color=Live)) +
  geom_point() +
  facet_wrap(~ Measure, ncol=2) + ggtitle("Phospho vs. Survival") + theme_bw()

# # LM Main
# model = lm(dss.days~CDK4+ADAM10+FGA+PAEP+HMOX1+CTNND1+PIK3CB, data=prot.clinical)
# xxx = summary(model)
# print(xxx)
# write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/proteome_lm.csv")
# 
# model = lm(OS~ADAM10_M+ADAM10_S+PIK3CB_M+PIK3CB_S+PAEP_M+PAEP_S+FGA_M+FGA_S+CDK4_M+CDK4_S+CTNND1_M+CTNND1_S+HMOX1_M+HMOX1_S, data=IHC)
# xxx = summary(model)
# print(xxx)
# write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/IHC_lm.csv")
# 
# # Logistic 
# model = glm(survival_6mo~CDK4+ADAM10+FGA+PAEP+HMOX1+CTNND1+PIK3CB, data=prot.clinical, family="binomial")
# xxx = summary(model)
# print(xxx)
# write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/proteome_logistic.csv")
# 
# model = glm(Live~ADAM10_M+ADAM10_S+PIK3CB_M+PIK3CB_S+PAEP_M+PAEP_S+FGA_M+FGA_S+CDK4_M+CDK4_S+CTNND1_M+CTNND1_S+HMOX1_M+HMOX1_S, data=IHC, family="binomial")
# xxx = summary(model)
# print(xxx)
# write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/IHC_logistic.csv")


# Cox Regression
IHC.cox = IHC
IHC.cox$status = abs(IHC.cox$Live-1)+1
IHC.cox$time = IHC.cox$OS
prot.clinical.cox = prot.clinical
prot.clinical.cox$status = prot.clinical.cox$dss.events+1
prot.clinical.cox$time = prot.clinical.cox$dss.days
phos.clinical.cox = phos.clinical
phos.clinical.cox$status = phos.clinical.cox$dss.events+1
phos.clinical.cox$time = phos.clinical.cox$dss.days

res.cox <- coxph(Surv(time, status) ~ ADAM10_M+ADAM10_S+PIK3CB_M+PIK3CB_S+PAEP_M+PAEP_S+FGA_M+FGA_S+CDK4_M+CDK4_S+CTNND1_M+CTNND1_S+HMOX1_M+HMOX1_S+NBP1_M+NBP1_S+DDX11_M, data =  IHC.cox)
aaa = summary(res.cox)
print(aaa)
write.csv(data.frame(aaa['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/IHC_cox.csv")
        
   
res.cox <- coxph(Surv(time, status) ~ CDK4+ADAM10+FGA+PAEP+HMOX1+CTNND1+PIK3CB+SCAI+DDX11, data=prot.clinical.cox)
xxx = summary(res.cox)
print(xxx)
write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/proteome_cox.csv")


res.cox <- coxph(Surv(time, status) ~ ADAM10+FGA+HMOX1+CTNND1, data=phos.clinical.cox)
ppp = summary(res.cox)
print(ppp)
write.csv(data.frame(ppp['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/phospho_cox.csv")

# prot
# Melanoma cells
ihc.coeff = data.frame(aaa['coefficients'])
ihc.coeff = ihc.coeff[grepl('_M', row.names(ihc.coeff)),]
ihc.coeff['gene'] = gsub('_M', '', row.names(ihc.coeff))
ihc.coeff$gene = gsub('NBP1', 'SCAI', ihc.coeff$gene)
row.names(ihc.coeff) = NULL
ihc.coeff = ihc.coeff[, c('gene', 'coefficients.z')]
colnames(ihc.coeff) = c('gene', 'IHC coef.')


prot.coeff = data.frame(xxx['coefficients'])
prot.coeff['gene'] = row.names(prot.coeff)
row.names(prot.coeff) = NULL
prot.coeff = prot.coeff[, c('gene', 'coefficients.z')]
colnames(prot.coeff) = c('gene', 'Proteome coef.')


coef = inner_join(ihc.coeff, prot.coeff, by="gene")
coef$Direction = coef$`IHC coef.`*coef$`Proteome coef.`>0
coef$Direction = gsub(TRUE, 'Match', coef$Direction)
coef$Direction = gsub(FALSE, 'Not Match', coef$Direction)
coef[round(coef$`IHC coef.`, digits=1) == 0 | round(coef$`Proteome coef.`, digits=1) == 0, 'Direction'] = "Undecided"
  

ggplot(coef, aes(x=`IHC coef.`, y=`Proteome coef.`, label=gene)) +
  geom_point(color = dplyr::case_when(coef$Direction == 'Match' ~ "Red", 
                                      coef$Direction == 'Undecided' ~ "Purple", 
                                      coef$Direction == 'Not Match' ~ "Blue"), 
             size=3)+
  geom_label_repel(aes(label = gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + xlim(-4,4) + ylim(-4,4) + theme_bw()

# Stromal cells
ihc.coeff = data.frame(aaa['coefficients'])
ihc.coeff = ihc.coeff[grepl('_S', row.names(ihc.coeff)),]
ihc.coeff['gene'] = gsub('_S', '', row.names(ihc.coeff))
ihc.coeff$gene = gsub('NBP1', 'SCAI', ihc.coeff$gene)
row.names(ihc.coeff) = NULL
ihc.coeff = ihc.coeff[, c('gene', 'coefficients.z')]
colnames(ihc.coeff) = c('gene', 'IHC coef.')


prot.coeff = data.frame(xxx['coefficients'])
prot.coeff['gene'] = row.names(prot.coeff)
row.names(prot.coeff) = NULL
prot.coeff = prot.coeff[, c('gene', 'coefficients.z')]
colnames(prot.coeff) = c('gene', 'Proteome coef.')

coef = inner_join(ihc.coeff, prot.coeff, by="gene")
coef$Direction = coef$`IHC coef.`*coef$`Proteome coef.`>0
coef$Direction = gsub(TRUE, 'Match', coef$Direction)
coef$Direction = gsub(FALSE, 'Not Match', coef$Direction)
coef[round(coef$`IHC coef.`, digits=0) == 0 | round(coef$`Proteome coef.`, digits=0) == 0, 'Direction'] = "Undecided"

ggplot(coef, aes(x=`IHC coef.`, y=`Proteome coef.`, label=gene)) +
  geom_point(color = dplyr::case_when(coef$Direction == 'Match' ~ "Red", 
                                      coef$Direction == 'Undecided' ~ "Purple", 
                                      coef$Direction == 'Not Match' ~ "Blue"), 
             size=3)+
  geom_label_repel(aes(label = gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + xlim(-4,4) + ylim(-4,4)+ theme_bw()


# phos
# Melanoma cells
ihc.coeff = data.frame(aaa['coefficients'])
ihc.coeff = ihc.coeff[grepl('_M', row.names(ihc.coeff)),]
ihc.coeff['gene'] = gsub('_M', '', row.names(ihc.coeff))
ihc.coeff$gene = gsub('NBP1', 'SCAI', ihc.coeff$gene)
row.names(ihc.coeff) = NULL
ihc.coeff = ihc.coeff[, c('gene', 'coefficients.z')]
colnames(ihc.coeff) = c('gene', 'IHC coef.')

phos.coeff = data.frame(ppp['coefficients'])
phos.coeff['gene'] = row.names(phos.coeff)
row.names(phos.coeff) = NULL
phos.coeff = phos.coeff[, c('gene', 'coefficients.z')]
colnames(phos.coeff) = c('gene', 'phospho coef.')


coef = inner_join(ihc.coeff, phos.coeff, by="gene")
coef$Direction = coef$`IHC coef.`*coef$`phospho coef.`>0
coef$Direction = gsub(TRUE, 'Match', coef$Direction)
coef$Direction = gsub(FALSE, 'Not Match', coef$Direction)
coef[round(coef$`IHC coef.`, digits=0) == 0 | round(coef$`phospho coef.`, digits=0) == 0, 'Direction'] = "Undecided"

ggplot(coef, aes(x=`IHC coef.`, y=`phospho coef.`, label=gene)) +
  geom_point(color = dplyr::case_when(coef$Direction == 'Match' ~ "Red",
                                      coef$Direction == 'Undecided' ~ "Purple",
                                      coef$Direction == 'Not Match' ~ "Blue"), 
             size=3)+
  geom_label_repel(aes(label = gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + xlim(-4,4) + ylim(-4,4)+ theme_bw()

# Stromal cells
ihc.coeff = data.frame(aaa['coefficients'])
ihc.coeff = ihc.coeff[grepl('_S', row.names(ihc.coeff)),]
ihc.coeff['gene'] = gsub('_S', '', row.names(ihc.coeff))
ihc.coeff$gene = gsub('NBP1', 'SCAI', ihc.coeff$gene)
row.names(ihc.coeff) = NULL
ihc.coeff = ihc.coeff[, c('gene', 'coefficients.z')]
colnames(ihc.coeff) = c('gene', 'IHC coef.')


phos.coeff = data.frame(ppp['coefficients'])
phos.coeff['gene'] = row.names(phos.coeff)
row.names(phos.coeff) = NULL
phos.coeff = phos.coeff[, c('gene', 'coefficients.z')]
colnames(phos.coeff) = c('gene', 'phospho coef.')

coef = inner_join(ihc.coeff, phos.coeff, by="gene")
coef$Direction = coef$`IHC coef.`*coef$`phospho coef.`>0
coef$Direction = gsub(TRUE, 'Match', coef$Direction)
coef$Direction = gsub(FALSE, 'Not Match', coef$Direction)
coef[round(coef$`IHC coef.`, digits=0) == 0 | round(coef$`phospho coef.`, digits=0) == 0, 'Direction'] = "Undecided"

ggplot(coef, aes(x=`IHC coef.`, y=`phospho coef.`, label=gene)) +
  geom_point(color = dplyr::case_when(coef$Direction == 'Match' ~ "Red",
                                      coef$Direction == 'Undecided' ~ "Purple",
                                      coef$Direction == 'Not Match' ~ "Blue"), 
             size=3)+
  geom_label_repel(aes(label = gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + xlim(-4,4) + ylim(-4,4)+ theme_bw()





# IHC ROC and KM
DDX11 <- read_excel("~/Documents/Segundo_Melanoma/Data/DDX11_IHC.xlsx")
DDX11 = DDX11[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_DDX11_N")]
DDX11$SPSSmodell_DDX11_N = gsub(1, "high", DDX11$SPSSmodell_DDX11_N)
DDX11$SPSSmodell_DDX11_N = gsub(0, "low", DDX11$SPSSmodell_DDX11_N)
DDX11$Live = abs(DDX11$Live-1)
DDX11_list = list(DDX11, "SPSSmodell_DDX11_N", "SPSSmodell_DDX11_N", "DDX11_N", "DDX11_N")
NBP1 <- read_excel("~/Documents/Segundo_Melanoma/Data/NBP1_IHC.xls")
NBP1 = NBP1[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_SCAI_M", "SPSSmodell_SCAI_S")]
NBP1$SPSSmodell_SCAI_M = gsub(1, "high", NBP1$SPSSmodell_SCAI_M)
NBP1$SPSSmodell_SCAI_M = gsub(0, "low", NBP1$SPSSmodell_SCAI_M)
NBP1$SPSSmodell_SCAI_S = gsub(0, "high", NBP1$SPSSmodell_SCAI_S)
NBP1$SPSSmodell_SCAI_S = gsub(1, "low", NBP1$SPSSmodell_SCAI_S)
NBP1$Live = abs(NBP1$Live-1)
NBP_list = list(NBP1, "SPSSmodell_SCAI_M", "SPSSmodell_SCAI_S", "SCAI_M", "SCAI_S")
ADAM10 <- read_excel("~/Documents/Segundo_Melanoma/Data/ADAM10_IHC.xls")
ADAM10 = ADAM10[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_ADAM10_M", "SPSSmodell_ADAM10_S")]
ADAM10$SPSSmodell_ADAM10_M = gsub(1, "high", ADAM10$SPSSmodell_ADAM10_M)
ADAM10$SPSSmodell_ADAM10_M = gsub(0, "low", ADAM10$SPSSmodell_ADAM10_M)
ADAM10$SPSSmodell_ADAM10_S = gsub(0, "high", ADAM10$SPSSmodell_ADAM10_S)
ADAM10$SPSSmodell_ADAM10_S = gsub(1, "low", ADAM10$SPSSmodell_ADAM10_S)
ADAM10$Live = abs(ADAM10$Live-1)
ADAM10_list = list(ADAM10, "SPSSmodell_ADAM10_M", "SPSSmodell_ADAM10_S", "ADAM10_M", "ADAM10_S")
PIK3CB <- read_excel("~/Documents/Segundo_Melanoma/Data/PIK3CB_IHC.xls")
PIK3CB = PIK3CB[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_PIK3cB_M", "SPSSmodell_PIK3cB_S")]
PIK3CB$SPSSmodell_PIK3cB_M = gsub(1, "high", PIK3CB$SPSSmodell_PIK3cB_M)
PIK3CB$SPSSmodell_PIK3cB_M = gsub(0, "low", PIK3CB$SPSSmodell_PIK3cB_M)
PIK3CB$SPSSmodell_PIK3cB_S = gsub(1, "high", PIK3CB$SPSSmodell_PIK3cB_S)
PIK3CB$SPSSmodell_PIK3cB_S = gsub(0, "low", PIK3CB$SPSSmodell_PIK3cB_S)
PIK3CB$Live = abs(PIK3CB$Live-1)
PIK3CB_list = list(PIK3CB, "SPSSmodell_PIK3cB_M", "SPSSmodell_PIK3cB_S", "PIK3cB_M", "PIK3cB_S")
PAEP <- read_excel("~/Documents/Segundo_Melanoma/Data/PAEP_IHC.xls")
PAEP = PAEP[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_PAEP_M", "SPSSmodell_PAEP_S")]
PAEP$SPSSmodell_PAEP_M = gsub(0, "high", PAEP$SPSSmodell_PAEP_M)
PAEP$SPSSmodell_PAEP_M = gsub(1, "low", PAEP$SPSSmodell_PAEP_M)
PAEP$SPSSmodell_PAEP_S = gsub(0, "high", PAEP$SPSSmodell_PAEP_S)
PAEP$SPSSmodell_PAEP_S = gsub(1, "low", PAEP$SPSSmodell_PAEP_S)
PAEP$Live = abs(PAEP$Live-1)
PAEP_list = list(PAEP, "SPSSmodell_PAEP_M", "SPSSmodell_PAEP_S", "PAEP_M", "PAEP_S")
FGA <- read_excel("~/Documents/Segundo_Melanoma/Data/FGA_IHC.xls")
FGA = FGA[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_FGA_M", "SPSSmodell_FGA_S")]
FGA$SPSSmodell_FGA_M = gsub(1, "high", FGA$SPSSmodell_FGA_M)
FGA$SPSSmodell_FGA_M = gsub(0, "low", FGA$SPSSmodell_FGA_M)
FGA$SPSSmodell_FGA_S = gsub(1, "high", FGA$SPSSmodell_FGA_S)
FGA$SPSSmodell_FGA_S = gsub(0, "low", FGA$SPSSmodell_FGA_S)
FGA$Live = abs(FGA$Live-1)
FGA_list = list(FGA, "SPSSmodell_FGA_M", "SPSSmodell_FGA_S", "FGA_M", "FGA_S")
CDK4 <- read_excel("~/Documents/Segundo_Melanoma/Data/CDK4_IHC.xls")
CDK4 = CDK4[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_CDK4_M", "SPSSmodell_CDK4_S")]
CDK4$SPSSmodell_CDK4_M = gsub(0, "high", CDK4$SPSSmodell_CDK4_M)
CDK4$SPSSmodell_CDK4_M = gsub(1, "low", CDK4$SPSSmodell_CDK4_M)
CDK4$SPSSmodell_CDK4_S = gsub(0, "high", CDK4$SPSSmodell_CDK4_S)
CDK4$SPSSmodell_CDK4_S = gsub(1, "low", CDK4$SPSSmodell_CDK4_S)
CDK4$Live = abs(CDK4$Live-1)
CDK4_list = list(CDK4, "SPSSmodell_CDK4_M", "SPSSmodell_CDK4_S", "CDK4_M", "CDK4_S")
HMOX1 <- read_excel("~/Documents/Segundo_Melanoma/Data/HMOX1_IHC.xlsx")
HMOX1 = HMOX1[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_HMOX_M", "SPSSmodell_HMOX_S")]
HMOX1$SPSSmodell_HMOX_M = gsub(0, "high", HMOX1$SPSSmodell_HMOX_M)
HMOX1$SPSSmodell_HMOX_M = gsub(1, "low", HMOX1$SPSSmodell_HMOX_M)
HMOX1$SPSSmodell_HMOX_S = gsub(0, "high", HMOX1$SPSSmodell_HMOX_S)
HMOX1$SPSSmodell_HMOX_S = gsub(1, "low", HMOX1$SPSSmodell_HMOX_S)
HMOX1$Live = abs(HMOX1$Live-1)
HMOX1_list = list(HMOX1, "SPSSmodell_HMOX_M", "SPSSmodell_HMOX_S", "HMOX_M", "HMOX_S")
CTNND1 <- read_excel("~/Documents/Segundo_Melanoma/Data/CTNND1_IHC.xls")
CTNND1 = CTNND1[, c("DFS", "PFS", "OS", "Live", "SPSSmodell_CTTND1_M", "SPSSmodell_CTTND1_S")]
CTNND1$SPSSmodell_CTTND1_M = gsub(0, "high", CTNND1$SPSSmodell_CTTND1_M)
CTNND1$SPSSmodell_CTTND1_M = gsub(1, "low", CTNND1$SPSSmodell_CTTND1_M)
CTNND1$SPSSmodell_CTTND1_S = gsub(0, "high", CTNND1$SPSSmodell_CTTND1_S)
CTNND1$SPSSmodell_CTTND1_S = gsub(1, "low", CTNND1$SPSSmodell_CTTND1_S)
CTNND1$Live = abs(CTNND1$Live-1)
CTNND1_list = list(CTNND1, "SPSSmodell_CTTND1_M", "SPSSmodell_CTTND1_S", "CTTND1_M", "CTTND1_S")
full_list = list(DDX11_list, NBP_list, ADAM10_list, PIK3CB_list, PAEP_list, FGA_list, CDK4_list, HMOX1_list, CTNND1_list)

for (i in full_list){
  for (s in c("OS", "DFS", "PFS")){
    km_fit <- survfit(Surv(get(s), Live) ~ get(i[[2]][1]), data=i[[1]])
    
    pv <- survminer::surv_pvalue(km_fit, method = "log-rank", test.for.trend = F, combine = F, data = i[[1]])  
    
    png(filename = paste("Results/IHC/km_", i[[4]][1], "_", s, ".png",sep =""), width = 500, height = 400, units = "px", pointsize = 16, bg = "white")
    km.plot <- autoplot(km_fit, conf.int = F, censor = F, surv.size = 2.5,  main = paste("p=", pv$pval, sep=""),
                        xlab = "",
                        ylab = "") +
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size=12, face="bold"), axis.text.y = element_text(size=12, face="bold"))
    
    print(km.plot)
    dev.off()
    
    km_fit <- survfit(Surv(get(s), Live) ~ get(i[[3]][1]), data=i[[1]])
    
    pv <- survminer::surv_pvalue(km_fit, method = "log-rank", test.for.trend = F, combine = F, data = i[[1]])  
    pv$pval
    
    png(filename = paste("Results/IHC/km_", i[[5]][1], "_", s, ".png",sep =""), width = 600, height = 400, units = "px", pointsize = 16, bg = "white")
    km.plot <- autoplot(km_fit, conf.int = F, censor = F, surv.size = 2.5,  main = paste("p=", pv$pval, sep=""),
                        xlab = "",
                        ylab = "") +
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
             axis.text.x = element_text(size=12, face="bold"), axis.text.y = element_text(size=12, face="bold"))
    
    print(km.plot)
    dev.off()
  }
}





