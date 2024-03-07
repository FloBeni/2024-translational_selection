# Generate Supplementary Figure 6, 7, 8 & 9
source("figure/figure_main_generator/library_path.R")

set_color_class = c("#E31A1C","#33A02C")
names(set_color_class) = c( "PO>nPO","nPO>PO")
set_shape_class = c(CDS=21,intron=24)
set_linetype_class = c(CDS="solid",intron="dashed")
set_alpha_class = c(CDS=1,intron=0.5)


###

data_11 = read.delim("data/data11_supp.tab")

dt_graph = data_11[ data_11$group == "POCs" ,]


# Pannel B

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_codon),big.mark=",",scientific=T)," SNPs PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_codon),big.mark=",",scientific=T)," SNPs nPO>PO",sep=""))

pB = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")

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
  axis.title.x = element_text(color="black", size=0,family="economica"),
  axis.title.y = element_text(color="black", size=30, family="economica"),
  axis.text.y =  element_text(color="black", size=30, family="economica"),
  axis.text.x =  element_text(color="black", size=0, family="economica"),
  title =  element_text(color="black", size=0, family="economica"),
  legend.text =  element_text(color="black", size=25, family="economica")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("SNP rate within CDSs")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F)+ ylim(0.01,0.055)+
  theme(legend.position = c(0.76, 0.9),
        legend.background = element_rect(fill="NA"),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.5, 'cm')
  ) + guides(fill= guide_legend(override.aes = list(size=7,pch=21), byrow = TRUE,order = 2))
pB



jpeg(paste(path_pannel,"p6B.jpg",sep=""), width = 6000/1, height = 2800/1,res=600/1)
print(pB)
dev.off()


# Pannel C

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_intron),big.mark=",",scientific=T)," SNPs PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_intron),big.mark=",",scientific=T)," SNPs nPO>PO",sep=""))

pC = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")

pC = pC +
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

pC = pC + theme_bw() + theme(
  axis.title.x = element_text(color="black", size=30,family="economica"),
  axis.title.y = element_text(color="black", size=30, family="economica"),
  axis.text.y =  element_text(color="black", size=30, family="economica"),
  axis.text.x =  element_text(color="black", size=27, family="economica"),
  title =  element_text(color="black", size=0, family="economica"),
  legend.text =  element_text(color="black", size=25, family="economica")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("SNP rate within introns")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F) + ylim(0.01,0.055) +
  theme(legend.position = c(0.76, 0.9),
        legend.background = element_rect(fill="NA"),
        legend.spacing.x = unit(.5, 'cm'),
        legend.spacing.y = unit(.5, 'cm')
  ) +  guides(fill= guide_legend(override.aes = list(size=7,pch=21),byrow = TRUE,order = 2))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
pC

jpeg(paste(path_pannel,"p6C.jpg",sep=""), width = 6000/1, height = 3000/1,res=600/1)
print(pC)
dev.off()


###

# data_10 = read.delim("data/data10_simulans_supp.tab")
data_10 = read.delim("data/data10_supp.tab")

dt_graph = data_10[ data_10$group == "POCs" ,]


# Pannel E

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_codon),big.mark=",",scientific=T)," substitutions PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_codon),big.mark=",",scientific=T)," substitutions nPO>PO",sep=""))

pE = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")

pE = pE +
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

pE = pE + theme_bw() + theme(
  axis.title.x = element_text(color="black", size=0,family="economica"),
  axis.title.y = element_text(color="black", size=30, family="economica"),
  axis.text.y =  element_text(color="black", size=30, family="economica"),
  axis.text.x =  element_text(color="black", size=0, family="economica"),
  title =  element_text(color="black", size=0, family="economica"),
  legend.text =  element_text(color="black", size=25, family="economica")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("Substitution rate within CDSs")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F) + ylim(0.005,0.05)+
  theme(legend.position = c(0.7, 0.9),
        legend.background = element_rect(fill="NA"),
        legend.spacing.x = unit(.5, 'cm'),
        legend.spacing.y = unit(.5, 'cm')
  ) +  guides(fill= guide_legend(override.aes = list(size=7,pch=21),byrow = TRUE,order = 2))
pE

jpeg(paste(path_pannel,"p6E.jpg",sep=""), width = 6000/1, height = 2800/1,res=600/1)
print(pE)
dev.off()


# Pannel F

labels_name = c( "PO>nPO" = paste( format(sum(dt_graph$sum_subst_density_optimal_to_nonoptimal_intron),big.mark=",",scientific=T)," substitutions PO>nPO",sep=""),
                 "nPO>PO"=paste( format(sum(dt_graph$sum_subst_density_nonoptimal_to_optimal_intron),big.mark=",",scientific=T)," substitutions nPO>PO",sep=""))

pF = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("POCs")

pF = pF +
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

pF = pF + theme_bw() + theme(
  axis.title.x = element_text(color="black", size=30,family="economica"),
  axis.title.y = element_text(color="black", size=30, family="economica"),
  axis.text.y =  element_text(color="black", size=30, family="economica"),
  axis.text.x =  element_text(color="black", size=27, family="economica"),
  title =  element_text(color="black", size=0, family="economica"),
  legend.text =  element_text(color="black", size=25, family="economica")
) + xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(limits = c(0.01,1000),
                breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
  ylab(paste("Substitution rate within introns")) + scale_fill_manual("",values=set_color_class,label=labels_name) + 
  scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
  scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
  annotation_logticks(sides="b")+
  guides(color = F, size = F , linetype=F, pch=F, alpha=F)  + ylim(0.005,0.05) +
  theme(legend.position = c(0.7, 0.9),
        legend.background = element_rect(fill="NA"),
        legend.spacing.x = unit(.5, 'cm'),
        legend.spacing.y = unit(.5, 'cm')
  ) +  guides(fill= guide_legend(override.aes = list(size=7,pch=21),byrow = TRUE,order = 2))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
pF

jpeg(paste(path_pannel,"p6F.jpg",sep=""), width = 6000/1, height = 3000/1,res=600/1)
print(pF)
dev.off()


# Figure 6

imgA = load.image(paste(path_require,"poly_subst.png",sep="") )
imgB = load.image(paste(path_pannel,"p6B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p6C.jpg",sep="") )
# imgD = load.image(paste(path_require,"substitutions.png",sep="") )
imgE = load.image(paste(path_pannel,"p6E.jpg",sep="") )
imgF = load.image(paste(path_pannel,"p6F.jpg",sep="") )
fly<-readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure6.pdf",sep=""), width=8, height=7)
  m = matrix(rep(NA,10*120), nrow=120)
  # for(i in 1:6){
  #   m[i,]=c(rep(1,10))
  # }
  
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
  par(mar=c(0, 1, 1.5, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  par(mar=c(0.5, 1, 0, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  
  par(mar=c(0, 1, 1.5, 0))
  plot(imgE, axes=FALSE)
  mtext("D",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  par(mar=c(0.5, 1, 0, 0))
  plot(imgF, axes=FALSE)
  mtext("E",at=-100,adj=-0, side=2, line=1, font=2, cex=1.6,las=2)
  
  dev.off()
}
