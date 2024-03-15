# Generate Figure 5
source("figure/figure_main_generator/library_path.R")


# Pannel 5 A

data6 = read.delim("data/data6_supp.tab")
data6$categorie = factor(data6$categorie,levels = rev( unique(data6$categorie))) 
dt_graph = data6[data6$species == "Homo_sapiens" & data6$type_aa == "POCs",]
# dt_graph = data6[data6$species == "Drosophila_melanogaster",]

pA = ggplot( dt_graph ,
              aes(y=freq,fill=categorie))  +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=set_color[c(4,3,5,6)]) +
  scale_shape_manual(values=c(24,22,21,23,25,20))+
  xlab("Constrain") + ylab("POC frequency") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=0, family="economica"),
    title =  element_text(color="black", size=18, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.42, face= "italic", size=18, family="economica"),
    plot.caption.position =  "plot"
  ) + scale_x_log10() + ggtitle(paste("N = " ,dt_graph$nb_genes[1]," BUSCO genes",sep="")) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  ) + theme(legend.position='none') + xlab("") + coord_cartesian(ylim=c(0.2,0.8))
pA



jpeg(paste(path_pannel,"p5A.jpg",sep=""),  width = 6200/2,  8000/2,res=900)
print(pA)
dev.off()


# Pannel 5 B

dt_graph = data6[data6$species == "Caenorhabditis_elegans" & data6$type_aa == "POCs",]


pB = ggplot( dt_graph ,
              aes(y=freq,fill=categorie))  +
  # geom_line(data=dt_graph[dt_graph$categorie %in% c("Highly constrained","Unconstrained"),],aes(group=busco_id),alpha=0.2,col="black",linetype="dashed")+
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual("Binning of codons per quartile\nof amino-acid constraint",values=set_color[c(4,3,5,6)],labels=c( "Unconstrained" = "25% least constrained" ,
                                                                        "Slighlty constrained" = "25-50%" ,
                                                                        "Moderately constrained" = "50-75%"  ,
                                                                        "Highly constrained" = "25% most constrained")) +
  scale_shape_manual(values=c(24,22,21,23,25,20)) +
  xlab("Constrain") + ylab("POC frequency") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=0, family="economica"),
    title =  element_text(color="black", size=18, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.42, face= "italic", size=18, family="economica"),
    plot.caption.position =  "plot"
  )  + ggtitle(paste("N = " ,dt_graph$nb_genes[1]," BUSCO genes",sep="")) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  )+ xlab("") + coord_cartesian(ylim=c(0.2,0.8)) +   
  theme(legend.spacing.y = unit(.4, 'cm'))   + 
  guides(fill = guide_legend(byrow = TRUE)) +ylab("")
pB

jpeg(paste(path_pannel,"p5B.jpg",sep=""),  width = 11000/2,  8000/2,res=900)
print(pB)
dev.off()


# Pannel 5 C

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]

pC = ggplot(data1,aes(y=constraint_overused_POCs,x=clade_group,fill=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=27, family="economica"),
    axis.text.y =  element_text(color="black", size=22, family="economica"),
    axis.text.x =  element_text(color="black",vjust=1,hjust=1, size=22,angle = 30, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + scale_fill_manual("Clades",values=Clade_color)+ 
  ylab(substitute(paste(Delta," POC"^"cons"))) + xlab("") + theme(legend.position='none')
# ylab("Difference in POC proportion between\nthe 25% highest and the 25% lowest constrained sites") +
pC


jpeg(paste(path_pannel,"p5C.jpg",sep=""), width = 5500/1, height = 3000/1,res=560/1)
print(pC)
dev.off()


# Figure 5

imgA = load.image(paste(path_pannel,"p5A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p5B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p5C.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans<-readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure5.pdf",sep=""), width=5.7, height=5)
  
  m = matrix(rep(NA,100*100), nrow=100)
  
  for(i in 1:50){
    m[i,]=c(rep(1,40),rep(2,60))
  }
  for(i in 50:100){
    m[i,]=c(rep(3,100))
  }
  layout(m)
  m
  
  par(mar=c(0, 0, 1, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=40,adj=-1, side=2, line=1, font=2, cex=1.3,las=2)
  xhuman=710
  yhuman=-460
  rasterImage(human,xleft=0+xhuman, ybottom=450/.8-yhuman, xright=190/.8+xhuman, ytop=0-yhuman)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=40,adj=-1, side=2, line=1, font=2, cex=1.3,las=2)
  xcel=720
  ycel=-480
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom = 350/1.6-ycel, xright = 1000/1.6+xcel, ytop=0-ycel)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgC, axes=FALSE)
  mtext("C",at=50,adj=-1, side=2, line=1, font=2, cex=1.3,las=2)
  dev.off()
}

