# Generate Supplementary Figure 4
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 4 A

data2_supp = read.delim(paste("data/data2_supp.tab",sep="") , header=T )
data2_supp$clade_group = GTDrift_list_species[data2_supp$species,]$clade_group

data2_supp$from = str_replace_all(data2_supp$from , "from_5prime","From start codons")
data2_supp$from = str_replace_all(data2_supp$from , "from_3prime","From stop codons")
data2_supp$from = factor(data2_supp$from,levels = c("From start codons" , "From stop codons"))


data2_supp[data2_supp$from == "From start codons",]$group = log10(data2_supp[data2_supp$from == "From start codons",]$group)
data2_supp[data2_supp$from == "From stop codons",]$group = -log10(data2_supp[data2_supp$from == "From stop codons",]$group)

species = "Homo_sapiens"
dt_graph = data2_supp[data2_supp$species == species,]

dt_graph$median_length = round(dt_graph$median_length / 1000 , 1)

dt_graph$median_length = factor(dt_graph$median_length,levels = rev(unique(dt_graph$median_length)))

pA = ggplot(dt_graph[dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,],aes(x = group , y= gci_count/posi_sites,fill=median_length,col=median_length)) +
  geom_line(lwd = 1.2,aes(x = group)) +
  # geom_point(pch=21,size=2,aes(x = group) )+
  facet_wrap(~from,scales = "free_x") +  theme_bw()  + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=26, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.30, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  ) + scale_fill_manual("Gene length (kb)",values = c("#1F78B4" , "#33A02C" , "#E31A1C" , "#FF7F00" ))+ scale_color_manual("Gene length (kb)",values = c("#1F78B4" , "#33A02C" , "#E31A1C" , "#FF7F00" )) + 
  scale_x_continuous(breaks=c(-1,-2,-3,-4,-5,1,2,3,4,5),labels = c(0.01,0.1,1,10,100,0.01,0.1,1,10,100)) +
  guides(fill = guide_legend(override.aes = list(lwd=3))) +
  # annotation_logticks(sides="b") +
  xlab("Distance from start codon or stop codon (kilobase, log scale)") + ylab ("GC introns") + ylim(c(0.3,0.7))
pA

jpeg(paste(path_pannel,"p4A_supp.jpg",sep=""), width = 5500/1, height = 2200/1,res=400/1)
print(pA)
dev.off()


# Supplementary Pannel 4 B

species = "Drosophila_melanogaster"
dt_graph = data2_supp[data2_supp$species == species,]

dt_graph$median_length = round(dt_graph$median_length / 1000,1)

dt_graph$median_length = factor(dt_graph$median_length,levels = rev(unique(dt_graph$median_length)))

pB = ggplot(dt_graph[dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,],aes(x = group,y= gci_count/posi_sites,fill=median_length,col=median_length)) +
  geom_line(data = dt_graph[ dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,] ,lwd = 1.2,aes(x = group,y= gci_count/posi_sites)) +
  # geom_point(data = dt_graph[ dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,] ,pch=21,size=2,aes(x = group,y= gci_count/posi_sites) )+
  facet_wrap(~from,scales = "free_x") +  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=26, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.30, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )  + scale_fill_manual("Gene length (kb)",values = c("#1F78B4" , "#33A02C" , "#E31A1C" , "#FF7F00" ))+ scale_color_manual("Gene length (kb)",values = c("#1F78B4" , "#33A02C" , "#E31A1C" , "#FF7F00" )) +
  scale_x_continuous(breaks=c(-1,-2,-3,-4,-5,1,2,3,4,5),labels = c(0.01,0.1,1,10,100,0.01,0.1,1,10,100)) +
  guides(fill = guide_legend(override.aes = list(lwd=3))) +
  # scale_x_log10(breaks=c(10,100,1000,10000),labels = c(10,100,"1kb","10kb")) +
  # annotation_logticks(sides="b") +
  xlab("Distance from start codon or stop codon (kilobase, log scale)") +  ylab ("GC introns") + ylim(c(0.3,0.7))
pB

jpeg(paste(path_pannel,"p4B_supp.jpg",sep=""), width = 5500/1, height = 2200/1,res=400/1)
print(pB)
dev.off()


# Supplementary Figure 4

imgA = load.image(paste(path_pannel,"p4A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p4B_supp.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
fly<-readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure4_supp.pdf",sep=""), width=4, height=3.5)
  
  m = matrix(rep(c(1,2),1*2), nrow=2)
  
  layout(m)
  m
  
  par(mar=c(0, 1, 1, 0))
  plot(imgA, axes=FALSE)
  mtext("A", side=2,at=-100,adj=-.5, line=1, font=2, cex=1,las=2)
  xhuman=400
  yhuman=-1300
  rasterImage(human,xleft=0+xhuman, ybottom=450/1-yhuman, xright=190/1+xhuman, ytop=0-yhuman)
  
  par(mar=c(0, 1, 1, 0))
  plot(imgB, axes=FALSE)
  mtext("B", side=2,at=-100,adj=-.5, line=1, font=2, cex=1,las=2)
  xaxis=400/1.1
  yaxis=1900/1.1
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1100/4+xaxis, ytop=-900/4+yaxis)
  dev.off()
}
