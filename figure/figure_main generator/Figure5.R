source("figure/figure_main generator/library_path.R")

############## Pannel 5 A
data5 = read.delim("data/data5_cons.tab")
data5 = data5[data5$type_aa == "Wb_WC_notambiguous",]
data5$categorie = factor(data5$categorie,levels = rev( unique(data5$categorie))) 
dt_graph = data5[data5$species == "Homo_sapiens",]
# dt_graph = data5[data5$species == "Drosophila_melanogaster",]

p5A = ggplot( dt_graph ,
              aes(y=freq,fill=categorie))  +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=set_color[c(4,3,5,6)]) +
  scale_shape_manual(values=c(24,22,21,23,25,20))+
  xlab("Constrain") + ylab("Frequency optimal codons") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=0, family="economica"),
    title =  element_text(color="black", size=22, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.42, face= "italic", size=18, family="economica"),
    plot.caption.position =  "plot"
  )  +
  scale_x_log10() + ggtitle(paste("Number of gene =" ,dt_graph$nb_genes[1])) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  ) + theme(legend.position='none') + xlab("") + coord_cartesian(ylim=c(0.2,0.8))
p5A



jpeg(paste(path_pannel,"p5A.jpg",sep=""),  width = 7000/2,  8000/2,res=1500/2)
print(p5A)
dev.off()



############## Pannel 5 B
dt_graph = data5[data5$species == "Caenorhabditis_elegans",]

p5B = ggplot( dt_graph ,
              aes(y=freq,fill=categorie))  +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual("Categories",values=set_color[c(4,3,5,6)]) +
  scale_shape_manual(values=c(24,22,21,23,25,20))+
  xlab("Constrain") + ylab("Frequency optimal codons") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=0, family="economica"),
    title =  element_text(color="black", size=22, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.42, face= "italic", size=18, family="economica"),
    plot.caption.position =  "plot"
  )  +
  scale_x_log10() + ggtitle(paste("Number of gene =" ,dt_graph$nb_genes[1])) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  )+ xlab("")+ coord_cartesian(ylim=c(0.2,0.8))
p5B

jpeg(paste(path_pannel,"p5B.jpg",sep=""),  width = 11000/2,  8000/2,res=1500/2)
print(p5B)
dev.off()


############## Pannel 5 C
data12 = read.delim("data/data12.tab")
data12 = data12[data12$species %in% clade_dt$species & !duplicated(data12$species),]
data12 = data12[ data12$nb_codon_not_decoded == 0  & data12$pval_aa_fpkm < 0.05 ,]

data6 = read.delim("data/data6.tab")
data6$clade_group = clade_dt[data6$species,]$clade_group

data6 = data6[data6$species %in% data12$species,]
data6 = data6[data6$type_aa == "Wb_WC_notambiguous",]
# data6 = data6[data6$nb_gene > 500,]

data6$ecart = data6$mean_highconstrain - data6$mean_unconstrain

p5C = ggplot(data6,aes(y=ecart*100,x=clade_group,fill=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0,angle = 50, family="economica"),
    title =  element_text(color="black", size=25, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    legend.position = "left"
  ) + scale_fill_manual("Clades",values=Clade_color)+ 
  ylab("Difference in proportion of optimal codons between\nthe 25% highest and the 25% lowest constrained sites") + xlab("")

p5C = ggMarginal(p5C, type="histogram",fill=set_color[1])
p5C


jpeg(paste(path_pannel,"p5C.jpg",sep=""), width = 7500/1, height = 3000/1,res=400/1)
print(p5C)
dev.off()







############## Figure 5

imgA = load.image(paste(path_pannel,"p5A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p5B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p5C.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans<-readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))
clade_png<-readPNG(paste(path_require,"clade.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure5.pdf",sep=""), width=7, height=5.5)
  
  m = matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:5){
    m[i,]=c(rep(1,5),rep(2,5))
  }
  for(i in 6:10){
    m[i,]=c(rep(3,10))
  }
  layout(m)
  m
  
  par(mar=c(1, 2, 1.7, 3))
  plot(imgA, axes=FALSE)
  mtext("A",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  xhuman=700
  yhuman=-500
  rasterImage(human,xleft=0+xhuman, ybottom=350/1-yhuman, xright=190/1+xhuman, ytop=0-yhuman)
  
  par(mar=c(1, 0, 1.7, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  xcel=700
  ycel=-400
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/2-ycel, xright=1000/2+xcel, ytop=0-ycel)
  
  par(mar=c(1, 0, 0, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  xmonkey=50
  ymonkey=500
  rasterImage(clade_png,xleft=0+xmonkey, ybottom=800/.45+ymonkey, xright=400/.45+xmonkey, ytop=ymonkey)
  dev.off()
}

