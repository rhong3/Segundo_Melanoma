# load ICAClusters.RData
# comparison with PCA results
library(pcaMethods)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(stringr)
library(cowplot)

cols = brewer.pal(6,name='Set1')
colside = droplevels(subtypes[colnames(proteome.data),'subtype'])
pr_pca = pca(t(proteome.data), nPcs=77,scale="uv", center=TRUE)
pca_scores = scores(pr_pca)
colnames(pca_scores)=1:77
colnames(pca_scores)=paste('PC',str_pad(as.character(colnames(pca_scores)),2,pad='0'),sep='_')
heatmap.2(t(pca_scores)[,colnames(proteome.data)[order(colside)]], trace ='none', scale='row',
          ColSideColors = cols[colside[order(colside)]],
          Rowv = F, 
          col=colorRampPalette(c("steelblue","white", "darkorange"))(64))
icP_pca=findCor2(t(pca_scores),clustering=1:77,clinical.data=clinical.data.num)
icP_pca=-log10(icP_pca)

icP_pca_plot = melt(icP_pca)
colnames(icP_pca_plot)=c('sig','clinical','value')
icP_pca_plot$sig=paste('PC',str_pad(as.character(icP_pca_plot$sig),2,pad='0'),sep='_')


p1=ggplot(icP_pca_plot,aes(x=sig,y=clinical))+
  theme_classic(base_size=12)+
  theme(axis.line=element_blank())+
  theme(axis.text.x = element_text(angle=90,vjust=0.5))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(low='white',high='dodgerblue')+
  geom_text(data=subset(icP_pca_plot,value>2),aes(label=round(value,digits=1)),size=2)+
  xlab('PCA-extracted signatures')+ylab('clinical features')+
  labs(fill='-log(P)')+
  ggtitle('Significance levels of correlation between clinical features and signature activity scores')

pr_ica = apply(ica50$A,2,function(x)tapply(x,as.factor(icCls$clustering),mean))
rownames(pr_ica)=paste('IC',str_pad(as.character(rownames(pr_ica)),2,pad='0'),sep='_')
heatmap.2(pr_ica[,colnames(proteome.data)[order(colside)]], trace ='none',scale ='none',
          ColSideColors = cols[colside[order(colside)]],
          Rowv=F,
          col=colorRampPalette(c("steelblue","white", "darkorange"))(64))
icP_ica=findCor2(pr_ica,clustering=1:77,clinical.data=clinical.data.num)
icP_ica=-log10(icP_ica)
icP_ica_plot = melt(icP_ica)
colnames(icP_ica_plot)=c('sig','clinical','value')
icP_ica_plot$sig=paste('IC',str_pad(as.character(icP_ica_plot$sig),2,pad='0'),sep='_')

p2=ggplot(icP_ica_plot,aes(x=sig,y=clinical))+
  theme_classic(base_size=12)+
  theme(axis.line=element_blank())+
  theme(axis.text.x = element_text(angle=90,vjust=0.5))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(low='white',high='dodgerblue')+
  geom_text(data=subset(icP_ica_plot,value>2),aes(label=round(value,digits=1)),size=2)+
  xlab('ICA-extracted signatures')+ylab('clinical features')+
  labs(fill='-log(P)')+
  ggtitle(' ')

plot_grid(p1,p2,ncol=1)

sum(apply(icP_pca,1,function(x)sum(x>3))>0)
sum(apply(icP_ica,1,function(x)sum(x>3))>0)
