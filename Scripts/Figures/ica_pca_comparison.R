# load ICAClusters.RData
# comparison with PCA results
library(pcaMethods)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(stringr)
library(cowplot)

ica_rdata = "Results/phospho/ICA/ICA_phospho_ICA.Rdata"
clinical_data = "Data/phospho/phospho_clinical.csv"
exp_prefix = "ICA_phospho"
feature = "phospho"

findCor2=function(mix,clustering,clinical.data,p_value=T){
  
  
  clsName=names(table(clustering))
  res=matrix(nrow=dim(mix)[1],ncol=dim(clinical.data)[2])
  colnames(res)=colnames(clinical.data)
  
  for (i in 1:dim(mix)[1]){
    
    for (j in 1:dim(clinical.data)[2]){
      
      #model.data=data.frame(y=clinical.data[,j],x=as.numeric(mix[i,]))
      #wprint(model.data)
      
      model.data=data.frame(y=as.numeric(clinical.data[,j]),x=as.numeric(mix[i,]))
      mod=lm(y~x,data=model.data)
      #print(colnames(clinical.data)[j])
      #print(summary(mod))
      if(p_value){
        res[i,j]=coef(summary(mod))[2,4] ## output can be p-value or coeffcients
      }
      else{
        res[i,j]=coef(summary(mod))[2,1]
      } 
    }
  }
  
  return(res)
}

clinical=read.table(file=clinical_data,
                    header=T,sep=',',row.names = 1)
#dim(clinical)

#ica_rdata='/media/lwk/lab/OV/retro/data/ov_pr_prosp_ICA.Rdata'
temp_space=new.env()
ica_vars_loaded=load(ica_rdata, temp_space)
#exp_prefix='ov_pr_prosp'


ica_vars=sub(".*_", "", ica_vars_loaded)

for (i in 1:length(ica_vars)){
  if(gsub(ica_vars[i],ica_vars_loaded[i],replacement = '')!=paste(exp_prefix,'_',sep='')){
    stop('Experiment name and ICA results data do not match. Double check the input!')
  }
  assign(ica_vars[i],get(ica_vars_loaded[i],envir=temp_space))
}


cols = brewer.pal(6,name='Set1')
# colside = droplevels(subtypes[colnames(dat),'subtype'])
pr_pca = pca(t(dat), nPcs=length(colnames(dat)),scale="uv", center=TRUE)
pca_scores = scores(pr_pca)
colnames(pca_scores)=1:length(colnames(dat))
colnames(pca_scores)=paste('PC',str_pad(as.character(colnames(pca_scores)),2,pad='0'),sep='_')
# heatmap.2(t(pca_scores)[,colnames(dat)], trace ='none', scale='row',
#           Rowv = F, 
#           col=colorRampPalette(c("steelblue","white", "darkorange"))(64))
# heatmap.2(t(pca_scores)[,colnames(dat)[order(colside)]], trace ='none', scale='row',
#           ColSideColors = cols[colside[order(colside)]],
#           Rowv = F, 
#           col=colorRampPalette(c("steelblue","white", "darkorange"))(64))
icP_pca=findCor2(t(pca_scores),clustering=1:length(colnames(dat)),clinical.data=clinical)
icP_pca=-log10(icP_pca)
row.names(icP_pca)=paste('PC',str_pad(as.character(seq(1:nrow(icP_pca))),2,pad='0'),sep='_')
row.names(icP_pca)[nchar(row.names(icP_pca))<6] = gsub("PC_", "PC_0", row.names(icP_pca)[nchar(row.names(icP_pca))<6])
icP_pca[is.na(icP_pca)] = 0
icP_pca_plot = melt(icP_pca)
colnames(icP_pca_plot)=c('sig','clinical','value')
icP_pca_plot=na.omit(icP_pca_plot)
icP_pca_plot=subset(icP_pca_plot,
                    sig%in%names(which(apply(icP_pca,1,max)>3))&clinical%in%colnames(icP_pca)[which(apply(icP_pca,2,max)>3)])

pr_ica = apply(ica$A,2,function(x)tapply(x,as.factor(iccls$clustering),mean))
# heatmap.2(pr_ica[,colnames(dat)[order(colside)]], trace ='none',scale ='none',
#           ColSideColors = cols[colside[order(colside)]],
#           Rowv=F,
#           col=colorRampPalette(c("steelblue","white", "darkorange"))(64))
icP_ica=findCor2(pr_ica,clustering=1:length(colnames(dat)),clinical.data=clinical)
icP_ica=-log10(icP_ica)
rownames(icP_ica)=paste('IC',str_pad(as.character(seq(1:nrow(icP_ica))),2,pad='0'),sep='_')
rownames(icP_ica)[nchar(rownames(icP_ica))<6] = gsub("IC_", "IC_0", rownames(icP_ica)[nchar(rownames(icP_ica))<6])
icP_ica[is.na(icP_ica)] = 0
icP_ica_plot = melt(icP_ica)
colnames(icP_ica_plot)=c('sig','clinical','value')
icP_ica_plot=na.omit(icP_ica_plot)
icP_ica_plot=subset(icP_ica_plot,
                   sig%in%names(which(apply(icP_ica,1,max)>3))&clinical%in%colnames(icP_ica)[which(apply(icP_ica,2,max)>3)])

rng = range(icP_ica_plot$value, icP_pca_plot$value)

p1=ggplot(icP_pca_plot,aes(x=sig,y=clinical))+
  theme_classic(base_size=20)+
  theme(axis.line=element_blank())+
  theme(text = element_text(size=20), axis.text.x = element_text(angle=90,vjust=0.5, size=20), axis.text.y=element_blank())+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(low='white',high='dodgerblue', midpoint=mean(rng)/1000, limits=c(floor(rng[1]), ceiling(rng[2])))+
  geom_text(data=subset(icP_pca_plot,value>3),aes(label=round(value,digits=1)),size=6)+
  xlab('PCA-extracted signatures (-log(P)>3)')+ylab('clinical and histological features')+
  labs(fill='-log(P)')+
  ggtitle(paste('Significance levels of correlation between clinical and histological features and', feature, 'data', sep=" "))

p2=ggplot(icP_ica_plot,aes(x=sig,y=clinical))+
  theme_classic(base_size=20)+
  theme(axis.line=element_blank())+
  theme(text = element_text(size=20), axis.text.x = element_text(angle=90,vjust=0.5, size=20), axis.text.y=element_blank(), legend.position="none")+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(low='white',high='dodgerblue', midpoint=mean(rng)/1000, limits=c(floor(rng[1]), ceiling(rng[2])))+
  geom_text(data=subset(icP_ica_plot,value>3),aes(label=round(value,digits=1)),size=6)+
  xlab('ICA-extracted signatures (-log(P)>3)')+ylab('clinical and histological features')+
  labs(fill='-log(P)')+
  ggtitle(' ')

png(filename=paste('~/documents/Segundo_Melanoma/Results/Figures/Figure3-S3/', feature, '_ICA_PCA.png', sep=''),
    width = 1200, height = 1000, units = "px", 
    pointsize = 20,
    bg = 'transparent')
plot_grid(p1,p2,ncol=1)
graphics.off()

