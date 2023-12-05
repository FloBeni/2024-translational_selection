# Generate Supplementary Figure 5
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 5 A

data3_supp = read.delim(paste("data/data9_supp.tab",sep="") , header=T )

data3_supp$clade_group = GTDrift_list_species[data3_supp$species,]$clade_group
data3_supp$from = str_replace_all(data3_supp$from , "from_5prime","From start codons")
data3_supp$from = str_replace_all(data3_supp$from , "from_3prime","From stop codons")
data3_supp$from = factor(data3_supp$from,levels = c("From start codons" , "From stop codons"))

data3_supp$label_facet = paste(data3_supp$clade_group, ", N=",  table(data3_supp[!duplicated(data3_supp$species),]$clade_group)[data3_supp$clade_group])
vector_label = unique(data3_supp$label_facet)
names(vector_label) = unique(data3_supp$clade_group)

data3_supp$label_facet = factor(data3_supp$label_facet , levels = vector_label[levels(data3_supp$clade_group)] )

data3_supp[data3_supp$from == "From start codons",]$group = log10(data3_supp[data3_supp$from == "From start codons",]$group)
data3_supp[data3_supp$from == "From stop codons",]$group = -log10(data3_supp[data3_supp$from == "From stop codons",]$group)

pA = ggplot(data3_supp[data3_supp$nb_genes >= 200,],aes(x = group,group=clade_group,y= gci_count/posi_sites  , label=nb_genes,fill=clade_group,col=clade_group)) +
  geom_line(lwd = 0.5,aes(group=species), alpha=0.7) + 
  facet_wrap(~label_facet + from , ncol=6 ,scale="free_x") +
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=15, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(family = "economica", size = 14)
  ) + scale_fill_manual(values = Clade_color)+ scale_color_manual(values = Clade_color) + 
  scale_x_continuous(breaks=c(-1,-2,-3,-4,1,2,3,4),labels = c(0.01,0.1,1,10,0.01,0.1,1,10)) +
  xlab("Distance from start codon or stop codon (kilobase, log scale)") + ylab ("GC introns")  + theme(legend.position='none') 
pA


jpeg(paste(path_pannel,"p5A_supp.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(pA)
dev.off()


# Supplementary Figure 5

imgA = load.image(paste(path_pannel,"p5A_supp.jpg",sep="") )

human<-readPNG(paste(path_require,"human.png",sep=""))

aves<-readPNG(paste(path_require,"aves.png",sep=""))
teleostei<-readPNG(paste(path_require,"teleostei.png",sep=""))
monkey<-readPNG(paste(path_require,"monkey.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))
bee<-readPNG(paste(path_require,"bee.png",sep=""))
tetrapod<-readPNG(paste(path_require,"tetrapod.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))
cnidaria = readPNG(paste(path_require,"cnidaria.png",sep=""))
insect = readPNG(paste(path_require,"insect.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure5_supp.pdf",sep=""), width=6.75, height=4)
  
  par(mar=c(0, 0, 0, 0))
  plot(imgA, axes=FALSE)
  # mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  
  xaxis=1000
  yaxis=2350
  rasterImage(monkey,xleft=0+xaxis, ybottom=0+yaxis, xright=900/5+xaxis, ytop=-900/5+yaxis)
  
  xaxis=2670
  yaxis=2350
  rasterImage(aves,xleft=0+xaxis, ybottom=0+yaxis, xright=600/3+xaxis, ytop=-750/3+yaxis)
  
  xaxis=4440
  yaxis=2380
  rasterImage(tetrapod,xleft=0+xaxis, ybottom=0+yaxis, xright=600/3.8+xaxis, ytop=-750/3.8+yaxis)
  
  xcel=900
  ycel=1350
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=0+ycel, xright=1000/4+xcel, ytop=-350/4+ycel)
  
  xaxis=4370
  yaxis=1350
  rasterImage(teleostei,xleft=0+xaxis, ybottom=0+yaxis, xright=1100/5+xaxis, ytop=-500/5+yaxis)
  
  xaxis=2670
  yaxis=450
  rasterImage(bee,xleft=0+xaxis, ybottom=0+yaxis, xright=900/5+xaxis, ytop=-700/5+yaxis)
  
  xaxis=970
  yaxis=450
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/8+xaxis, ytop=-900/8+yaxis)
  
  xaxis=4430
  yaxis=450
  rasterImage(insect,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/7+xaxis, ytop=-1100/7+yaxis)
  
  xcel=2640
  ycel=1410
  rasterImage(cnidaria,xleft=0+xcel, ybottom=0+ycel, xright=900/4+xcel, ytop=-750/4+ycel)
  dev.off()
}

