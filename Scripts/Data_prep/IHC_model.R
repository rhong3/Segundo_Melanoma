# IHC survival models
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library("survival")
library("survminer")


# IHC data filter
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
  select(c(1:3, 8:19, 4:7))

write.csv(IHC, "~/Documents/Segundo_Melanoma/Data/IHC_agg.csv", row.names = FALSE)
  

# Proteomics data filter
gene = c('PIK3CB', 'PAEP', 'FGA', 'CDK4', 'ADAM10', 'HMOX1', 'CTNND1')
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


# LM 
# prep
prot.clinical = read.csv("~/Documents/Segundo_Melanoma/Data/IHC_prot_clinical.csv")
IHC = read.csv("~/Documents/Segundo_Melanoma/Data/IHC_agg.csv")
colnames(IHC) = gsub("PIK3cB", "PIK3CB", colnames(IHC))
colnames(IHC) = gsub("HMOX", "HMOX1", colnames(IHC))
gene_means = colMeans(prot.clinical[,2:8])
gene_std = apply(prot.clinical[,2:8], 2, sd)
IHC.scale = IHC
center_factors = c(0.0302873, 0.0302873, 0.0896912, 0.0896912, -0.2300511, -0.2300511, -0.1936995, -0.1936995, 0.1471385, 0.1471385, -0.2320132, -0.2320132, 0.1701781, 0.1701781)
scale_factors = c(0.3197894, 0.3197894, 0.4426608, 0.4426608, 1.4249570, 1.4249570, 1.3063025, 1.3063025, 0.5791502, 0.5791502, 0.5169386, 0.5169386, 0.6680108, 0.6680108)
IHC.scale[, 2:15] = scale(IHC.scale[, 2:15], center=scale_factors, scale=FALSE)

# plot
IHC.plot = data.frame()
for (i in 2:15){
  IHC.plot.temp = data.frame(Values=IHC[,i], Measure=colnames(IHC)[i], Days=IHC$OS, Live=IHC$Live)
  IHC.plot = rbind.data.frame(IHC.plot, IHC.plot.temp)
}
IHC.plot$Live = as.factor(IHC.plot$Live)
ggplot(IHC.plot, aes(Days, Values, color=Live)) +
  geom_point() +
  facet_wrap(~ Measure, ncol=7) + ggtitle("IHC vs. Survival") + theme_bw()

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
for (i in 2:8){
  prot.clinical.temp = data.frame(Values=prot.clinical[,i], Measure=colnames(prot.clinical)[i], Days=prot.clinical$dss.days, Live=abs(prot.clinical$dss.events-1))
  prot.clinical.plot = rbind.data.frame(prot.clinical.plot, prot.clinical.temp)
}
prot.clinical.plot$Live = as.factor(prot.clinical.plot$Live)
ggplot(prot.clinical.plot, aes(Days, Values, color=Live)) +
  geom_point() +
  facet_wrap(~ Measure, ncol=4) + ggtitle("Proteome vs. Survival") + theme_bw()

# LM Main
model = lm(dss.days~CDK4+ADAM10+FGA+PAEP+HMOX1+CTNND1+PIK3CB, data=prot.clinical)
xxx = summary(model)
print(xxx)
write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/proteome_lm.csv")

model = lm(OS~ADAM10_M+ADAM10_S+PIK3CB_M+PIK3CB_S+PAEP_M+PAEP_S+FGA_M+FGA_S+CDK4_M+CDK4_S+CTNND1_M+CTNND1_S+HMOX1_M+HMOX1_S, data=IHC)
xxx = summary(model)
print(xxx)
write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/IHC_lm.csv")

# Logistic 
model = glm(survival_6mo~CDK4+ADAM10+FGA+PAEP+HMOX1+CTNND1+PIK3CB, data=prot.clinical, family="binomial")
xxx = summary(model)
print(xxx)
write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/proteome_logistic.csv")

model = glm(Live~ADAM10_M+ADAM10_S+PIK3CB_M+PIK3CB_S+PAEP_M+PAEP_S+FGA_M+FGA_S+CDK4_M+CDK4_S+CTNND1_M+CTNND1_S+HMOX1_M+HMOX1_S, data=IHC, family="binomial")
xxx = summary(model)
print(xxx)
write.csv(data.frame(xxx['coefficients']),  "~/Documents/Segundo_Melanoma/Results/IHC/IHC_logistic.csv")


# Cox Regression



