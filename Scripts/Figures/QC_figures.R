

#### Box plot for sequence coverage B1-B15

TMT11_seq_cov_B1B15 <- read.table("Data/others/Secuence_coverage_in_percent_B1-B15.txt", header = TRUE)

boxplot(Sequence_coverage ~ Batch, data = TMT11_seq_cov_B1B15)

Batches <- factor(TMT11_seq_cov_B1B15$Batch , levels=c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15")) #Give the boxes a specific order


boxplot(Sequence_coverage ~ Batches, data = TMT11_seq_cov_B1B15, col = "#a9afd3", ylab = "Sequence coverage per protein (%)", xlab = "TMT11-plex experiment", cex.lab=1.5, cex.axis=1.5)




#### Bland-Altman plot with ggplot2 TMT11 B1 ----

library(BlandAltmanLeh)
library(viridis)
library(ggplot2)
library(RColorBrewer)

TMT11_MM_b1_ <- read.table("Data/others/20200421_B1_log2_100vv_log2_Bland-Altman_plot.txt", header = TRUE)


bland.altman.plot(TMT11_MM_b1_$R1, TMT11_MM_b1_$R2, xlab="Means (TMT11 B1 R1-R2)", ylab="Difference", head= "MM TMT11 B1", silent=FALSE)

p <- ggplot(TMT11_MM_b1_, aes(x=Mean,y=Dif))+geom_point(shape=21,size=2,alpha=0.5,aes(fill=Dif))+geom_hline(yintercept = c( -3.363, -1.2969, 0.7695067), size=1,linetype=c("dashed", "solid", "dashed"),color=c("black", "red", "black"))+theme_bw()+scale_fill_viridis(option="plasma",begin=0.2,end=0.8)
p + ggtitle("TMT11 B1") + xlab("Mean (log2 intensity)") + ylab("Difference (R1-R2)")+theme(plot.title = element_text(size=22, face="bold", hjust = 0.5), axis.text=element_text(size=22), axis.title = element_text(size = 22))


## blue (Safe 16 SVG Hex3)	#0000FF	
#9be1fb; #a1cae6; #a9afd3; #a8a7c0; #f799d1; #f9b3af; #ffeca1; #f4b4ad; #fffbb2
