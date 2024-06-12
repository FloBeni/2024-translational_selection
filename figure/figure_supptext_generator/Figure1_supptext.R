# Generate Supplementary Texts Figure 1
source("figure/figure_supptext_generator/library_path.R")

set_color_class = c("#E31A1C","#33A02C")
names(set_color_class) = c( "PO>nPO","nPO>PO")
set_shape_class = c(CDS=21,intron=24)
set_linetype_class = c(CDS="solid",intron="dashed")
set_alpha_class = c(CDS=1,intron=0.5)
resolution = 4

# Pannel B

data_11 = read.delim("data/data11_supp.tab")
dt_graph = data_11

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_codon),big.mark=",",scientific=T)," SNPs PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_codon),big.mark=",",scientific=T)," SNPs nPO>PO",sep=""))


pB = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")+ 
  geom_text(data=data.frame(),label="CDS",aes(x=0.012,y=0.053), family="ubuntu condensed", size=12)

pB = pB +
  geom_line(aes(y=density_optimal_to_nonoptimal_codon,color="PO>nPO",linetype='CDS',alpha='CDS'),size=1) +
  geom_point(aes(y=density_optimal_to_nonoptimal_codon,fill="PO>nPO",pch='CDS',alpha='CDS'),size=3)+

  geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_codon,
                    ymax=confint_high_density_optimal_to_nonoptimal_codon,color="PO>nPO"),
                width=0.03,show.legend=FALSE) +
  
  geom_line(aes(y=density_nonoptimal_to_optimal_codon,color="nPO>PO",linetype='CDS',alpha='CDS'),size=1) +
  geom_point(aes(y=density_nonoptimal_to_optimal_codon,fill="nPO>PO",pch='CDS',alpha='CDS'),size=3)+

  geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_codon,
                    ymax=confint_high_density_nonoptimal_to_optimal_codon,color="nPO>PO"),
                width=0.03,show.legend=FALSE)

pB = pB + theme_bw() + theme(
  axis.title.x = element_text(color="black",vjust=0.5, size=0,family="ubuntu condensed"),
  axis.title.y = element_text(color="black",vjust=1.5, size=30, family="ubuntu condensed"),
  axis.text.y =  element_text(color="black", size=30, family="ubuntu condensed"),
  axis.text.x =  element_text(color="black", size=0, family="ubuntu condensed"),
  title =  element_text(color="black", size=0, family="ubuntu condensed"),
  legend.text =  element_text(color="black", size=25, family="ubuntu condensed")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("SNP rate")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F)+ ylim(0.01,0.055)+
  theme(legend.position = c(0.76, 0.87),
        legend.background = element_rect(fill="NA"),
        legend.key.spacing.x = unit(0.5, 'cm'),
        legend.key.spacing.y = unit(.5, 'cm')
  ) + guides(fill= guide_legend(override.aes = list(size=7,pch=21), byrow = TRUE,order = 2)) 
pB



jpeg(paste(path_pannel,"p1B_supptext.jpg",sep=""), width = 6000/1/resolution, height = 2800/1/resolution,res=600/1/resolution)
print(pB)
dev.off()


# Pannel C

data_10 = read.delim("data/data10_supp.tab")
dt_graph = data_10

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_codon),big.mark=",",scientific=T)," substitutions PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_codon),big.mark=",",scientific=T)," substitutions nPO>PO",sep=""))

pC = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")+ 
  geom_text(data=data.frame(),label="CDS",aes(x=0.012,y=0.048), family="ubuntu condensed", size=12)

pC = pC +
  geom_line(aes(y=density_optimal_to_nonoptimal_codon,color="PO>nPO",linetype='CDS',alpha='CDS'),size=1) +
  geom_point(aes(y=density_optimal_to_nonoptimal_codon,fill="PO>nPO",pch='CDS',alpha='CDS'),size=3)+
  
  geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_codon,
                    ymax=confint_high_density_optimal_to_nonoptimal_codon,color="PO>nPO"),
                width=0.03,show.legend=FALSE) +
  
  
  geom_line(aes(y=density_nonoptimal_to_optimal_codon,color="nPO>PO",linetype='CDS',alpha='CDS'),size=1) +
  geom_point(aes(y=density_nonoptimal_to_optimal_codon,fill="nPO>PO",pch='CDS',alpha='CDS'),size=3)+
  geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_codon,
                    ymax=confint_high_density_nonoptimal_to_optimal_codon,color="nPO>PO"),
                width=0.03,show.legend=FALSE)

pC = pC + theme_bw() + theme(
  axis.title.x = element_text(color="black",vjust=0, size=0,family="ubuntu condensed"),
  axis.title.y = element_text(color="black",vjust=1.5, size=30, family="ubuntu condensed"),
  axis.text.y =  element_text(color="black", size=30, family="ubuntu condensed"),
  axis.text.x =  element_text(color="black", size=0, family="ubuntu condensed"),
  title =  element_text(color="black", size=0, family="ubuntu condensed"),
  legend.text =  element_text(color="black", size=25, family="ubuntu condensed")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("Substitution rate")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F) + ylim(0.005,0.05)+
  theme(legend.position = c(0.7, 0.87),
        legend.background = element_rect(fill="NA"),
        legend.key.spacing.x = unit(.5, 'cm'),
        legend.key.spacing.y = unit(.5, 'cm')
  ) +  guides(fill= guide_legend(override.aes = list(size=7,pch=21),byrow = TRUE,order = 2))
pC

jpeg(paste(path_pannel,"p1C_supptext.jpg",sep=""), width = 6000/1/resolution, height = 2800/1/resolution,res=600/1/resolution)
print(pC)
dev.off()



# Pannel D

data_11 = read.delim("data/data11_supp.tab")
dt_graph = data_11

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_intron),big.mark=",",scientific=T)," SNPs PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_intron),big.mark=",",scientific=T)," SNPs nPO>PO",sep=""))

pD = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")+ 
  geom_text(data=data.frame(),label="Intron",aes(x=0.015,y=0.053), family="ubuntu condensed", size=12)

pD = pD +
  geom_line(aes(y=density_optimal_to_nonoptimal_intron,color="PO>nPO",linetype='intron',alpha='intron'),size=1) +
  geom_point(aes(y=density_optimal_to_nonoptimal_intron,fill="PO>nPO",pch='intron',alpha='intron'),size=3)+
  geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_intron,
                    ymax=confint_high_density_optimal_to_nonoptimal_intron,color="PO>nPO"),
                width=0.03,show.legend=FALSE) +
  
  geom_line(aes(y=density_nonoptimal_to_optimal_intron,color="nPO>PO",linetype='intron',alpha='intron'),size=1) +
  geom_point(aes(y=density_nonoptimal_to_optimal_intron,fill="nPO>PO",pch='intron',alpha='intron'),size=3)+
  geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_intron,
                    ymax=confint_high_density_nonoptimal_to_optimal_intron,color="nPO>PO"),
                width=0.03,show.legend=FALSE)

pD = pD + theme_bw() + theme(
  axis.title.x = element_text(color="black",vjust=-0, size=25,family="ubuntu condensed"),
  axis.title.y = element_text(color="black", size=30, family="ubuntu condensed"),
  axis.text.y =  element_text(color="black",vjust=1.5, size=30, family="ubuntu condensed"),
  axis.text.x =  element_text(color="black", size=27, family="ubuntu condensed"),
  title =  element_text(color="black", size=0, family="ubuntu condensed"),
  legend.text =  element_text(color="black", size=25, family="ubuntu condensed")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("SNP rate")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F) + ylim(0.01,0.055) +
  theme(legend.position = c(0.76, 0.87),
        legend.background = element_rect(fill="NA"),
        legend.key.spacing.x = unit(.5, 'cm'),
        legend.key.spacing.y = unit(.5, 'cm')
  ) +  guides(fill= guide_legend(override.aes = list(size=7,pch=21),byrow = TRUE,order = 2))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
pD

jpeg(paste(path_pannel,"p1D_supptext.jpg",sep=""), width = 6000/1/resolution, height = 3000/1/resolution,res=600/1/resolution)
print(pD)
dev.off()



# Pannel E

data_10 = read.delim("data/data10_supp.tab")
dt_graph = data_10

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_intron),big.mark=",",scientific=T)," substitutions PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_intron),big.mark=",",scientific=T)," substitutions nPO>PO",sep=""))

pE = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")+ 
  geom_text(data=data.frame(),label="Intron",aes(x=0.015,y=0.048), family="ubuntu condensed", size=12)

pE = pE +
  geom_line(aes(y=density_optimal_to_nonoptimal_intron,color="PO>nPO",linetype='intron',alpha='intron'),size=1) +
  geom_point(aes(y=density_optimal_to_nonoptimal_intron,fill="PO>nPO",pch='intron',alpha='intron'),size=3)+
  geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_intron,
                    ymax=confint_high_density_optimal_to_nonoptimal_intron,color="PO>nPO"),
                width=0.03,show.legend=FALSE) +
  
  geom_line(aes(y=density_nonoptimal_to_optimal_intron,color="nPO>PO",linetype='intron',alpha='intron'),size=1) +
  geom_point(aes(y=density_nonoptimal_to_optimal_intron,fill="nPO>PO",pch='intron',alpha='intron'),size=3)+
  geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_intron,
                    ymax=confint_high_density_nonoptimal_to_optimal_intron,color="nPO>PO"),
                width=0.03,show.legend=FALSE)

pE = pE + theme_bw() + theme(
  axis.title.x = element_text(color="black",vjust=0, size=25,family="ubuntu condensed"),
  axis.title.y = element_text(color="black",vjust=1.5, size=30, family="ubuntu condensed"),
  axis.text.y =  element_text(color="black", size=30, family="ubuntu condensed"),
  axis.text.x =  element_text(color="black", size=27, family="ubuntu condensed"),
  title =  element_text(color="black", size=0, family="ubuntu condensed"),
  legend.text =  element_text(color="black", size=25, family="ubuntu condensed")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("Substitution rate")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F)  + ylim(0.005,0.05) +
  theme(legend.position = c(0.7, 0.87),
        legend.background = element_rect(fill="NA"),
        legend.key.spacing.x = unit(.5, 'cm'),
        legend.key.spacing.y = unit(.5, 'cm')
  ) +  guides(fill= guide_legend(override.aes = list(size=7,pch=21),byrow = TRUE,order = 2))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
pE

jpeg(paste(path_pannel,"p1E_supptext.jpg",sep=""), width = 6000/1/resolution, height = 3000/1/resolution,res=600/1/resolution)
print(pE)
dev.off()


# Supplementary Texts Figure 1

imgA = load.image(paste(path_require,"poly_subst.png",sep="") )
imgB = load.image(paste(path_pannel,"p1B_supptext.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p1C_supptext.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p1D_supptext.jpg",sep="") )
imgE = load.image(paste(path_pannel,"p1E_supptext.jpg",sep="") )
fly<-readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure1_supptext.pdf",sep=""), width=8, height=7)
  m = matrix(rep(NA,10*120), nrow=120)
  
  for(i in 1:6){
    m[,i]=c(rep(1,50),rep(2,35),rep(3,35))
  }
  for(i in 6:10){
    m[,i]=c(rep(1,50),rep(4,35),rep(5,35))
  }
  m
  layout(m)
  
  par(mar=c(0, 0, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=100,adj=-4, side=2, line=1, font=2, cex=1.7,las=2)
  xaxis=1300/1.1
  yaxis=120/1.1
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1100/13+xaxis, ytop=-900/13+yaxis)
  par(mar=c(0, 1, 1.5, 1))
  plot(imgB, axes=FALSE)
  mtext("B",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  par(mar=c(0.5, 1, 0, 1))
  plot(imgD, axes=FALSE)
  mtext("D",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  
  par(mar=c(0, 1, 1.5, 1))
  plot(imgC, axes=FALSE)
  mtext("C",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  par(mar=c(0.5, 1, 0, 1))
  plot(imgE, axes=FALSE)
  mtext("E",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  
  dev.off()
}
