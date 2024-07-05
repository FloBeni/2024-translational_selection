# Generate Supplementary Texts Figure 3
source("figure/figure_supptext_generator/library_path.R")
resolution = 3

# Pannel A
data6 = read.delim("data/data6_supp.tab",comment.char = "#")
data6$categorie = factor(data6$categorie,levels = rev( unique(data6$categorie))) 
dt_graph = data6[data6$species == "Homo_sapiens" & data6$set == "POCs",]

pA = ggplot( dt_graph , aes(y=freq,fill=categorie))  +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=set_color[c(4,3,5,6)]) +
  scale_shape_manual(values=c(24,22,21,23,25,20))+
  xlab("Constrain") + ylab("POC frequency") + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=0.5, size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=25, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=20, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=0, family="ubuntu condensed"),
    title =  element_text(color="black", size=15, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.42, face= "italic", size=18, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  ) + ggtitle(paste("N = " ,dt_graph$nb_genes[1]," BUSCO genes",sep="")) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  ) + theme(legend.position='none') + xlab("") + coord_cartesian(ylim=c(0.2,0.8))
pA



jpeg(paste(path_pannel,"p3A_supptext.jpg",sep=""),  width = 6200/2/resolution,  8000/2/resolution,res=900/resolution)
print(pA)
dev.off()


# Pannel B

dt_graph = data6[data6$species == "Caenorhabditis_elegans" & data6$set == "POCs",]


pB = ggplot( dt_graph , aes(y=freq,fill=categorie))  +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual("Binning of codons per quartile\nof amino-acid constraint",values=set_color[c(4,3,5,6)],labels=c( "Unconstrained" = "25% least constrained" ,
                                                                        "Slighlty constrained" = "25-50%" ,
                                                                        "Moderately constrained" = "50-75%"  ,
                                                                        "Highly constrained" = "25% most constrained")) +
  scale_shape_manual(values=c(24,22,21,23,25,20)) +
  xlab("Constrain") + ylab("POC frequency") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=25, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=20, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=0, family="ubuntu condensed"),
    title =  element_text(color="black", size=15, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=16, family="ubuntu condensed"),
    legend.title =  element_text(color="black", size=17, family="ubuntu condensed"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.42, face= "italic", size=18, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )  + ggtitle(paste("N = " ,dt_graph$nb_genes[1]," BUSCO genes",sep="")) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  )+ xlab("") + coord_cartesian(ylim=c(0.2,0.8)) +   
  # theme(legend.spacing.y = unit(.4, 'cm'))   + 
  guides(fill = guide_legend(byrow = TRUE)) +ylab("")
pB

jpeg(paste(path_pannel,"p3B_supptext.jpg",sep=""),  width = 11500/2/resolution,  8000/2/resolution,res=900/resolution)
print(pB)
dev.off()


# Pannel C

data1 = read.delim("data/data1_supp.tab",comment.char = "#")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]

pC = ggplot(data1,aes(y=constraint_overused_POCs,x=clade_group,fill=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=2.5,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=27, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=22, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black",vjust=1,hjust=1, size=22,angle = 30, family="ubuntu condensed"),
    title =  element_text(color="black", size=15, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed")
  ) + scale_fill_manual("Clades",values=Clade_color)+ 
  ylab(substitute(paste(Delta," POC"^"cons"))) + xlab("") + theme(legend.position='none')
pC


jpeg(paste(path_pannel,"p3C_supptext.jpg",sep=""), width = 5500/1/resolution, height = 3000/1/resolution,res=560/1/resolution)
print(pC)
dev.off()


# Supplementary Texts Figure 3

imgA = load.image(paste(path_pannel,"p3A_supptext.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p3B_supptext.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p3C_supptext.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans<-readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure3_supptext.pdf",sep=""), width=5.7, height=5)
  
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
  xhuman=740/resolution
  yhuman=-460/resolution
  rasterImage(human,xleft=0+xhuman, ybottom=450/.8/resolution-yhuman, xright=190/.8/resolution+xhuman, ytop=0-yhuman)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=40,adj=-1, side=2, line=1, font=2, cex=1.3,las=2)
  xcel=740/resolution
  ycel=-480/resolution
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom = 350/1.6/resolution-ycel, xright = 1000/1.6/resolution+xcel, ytop=0-ycel)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgC, axes=FALSE)
  mtext("C",at=50,adj=-1, side=2, line=1, font=2, cex=1.3,las=2)
  dev.off()
}

