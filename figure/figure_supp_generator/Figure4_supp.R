source("figure/figure_supp generator/library_path.R")


############## Supplementary Pannel 4 A
arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label


data11 = read.delim(paste("data/data11.tab",sep="") , header=T )

data11$clade_group = clade_dt[data11$species,]$clade_group

data11$from = str_replace_all(data11$from , "from_5prime","From start")
data11$from = str_replace_all(data11$from , "from_3prime","From stop")
data11$from = factor(data11$from,levels = c("From start" , "From stop"))


data11[data11$from == "From start",]$group = log10(data11[data11$from == "From start",]$group)
data11[data11$from == "From stop",]$group = -log10(data11[data11$from == "From stop",]$group)

species = "Homo_sapiens"
dt_graph = data11[data11$species == species,]

dt_graph$median_length = round(dt_graph$median_length / 1000 , 1)

dt_graph$median_length = factor(dt_graph$median_length,levels = rev(unique(dt_graph$median_length)))

p4A = ggplot(dt_graph[dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,],aes(x = group,y= gci_count/posi_sites,fill=median_length,col=median_length)) +
  geom_line(data = dt_graph[ dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,] ,lwd = 1.2,aes(x = group,y= gci_count/posi_sites)) +
  geom_point(data = dt_graph[ dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,] ,pch=21,size=2,aes(x = group,y= gci_count/posi_sites) )+
  facet_wrap(~from,scales = "free_x") +  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=25, family="economica"),
    legend.text =  element_text(color="black", size=25, family="economica"),
    strip.text = element_text(family = "economica", size = 14)
  ) + scale_fill_manual("Gene length (kb)",values = set_color)+ scale_color_manual("Gene length (kb)",values = set_color) + ggtitle(species) +
  scale_x_continuous(breaks=c(-1,-2,-3,-4,-5,1,2,3,4,5),labels = c(0.01,0.1,1,10,100,0.01,0.1,1,10,100)) +
  # scale_x_log10(breaks=c(10,100,1000,10000),labels = c(10,100,"1kb","10kb")) +
  # annotation_logticks(sides="b") +
  xlab("Genome distance from start codon or stop codon (log10 scale)") + ylab ("GCi") + ylim(c(0.3,0.7))
print( p4A )

jpeg(paste(path_pannel,"p4A_supp.jpg",sep=""), width = 5500/1, height = 2500/1,res=400/1)
print(p4A)
dev.off()

############## Supplementary Pannel 4 B
species = "Drosophila_melanogaster"
dt_graph = data11[data11$species == species,]

dt_graph$median_length = round(dt_graph$median_length / 1000,1)

dt_graph$median_length = factor(dt_graph$median_length,levels = rev(unique(dt_graph$median_length)))

p4B = ggplot(dt_graph[dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,],aes(x = group,y= gci_count/posi_sites,fill=median_length,col=median_length)) +
  geom_line(data = dt_graph[ dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,] ,lwd = 1.2,aes(x = group,y= gci_count/posi_sites)) +
  geom_point(data = dt_graph[ dt_graph$nb_genes >= max(dt_graph$nb_genes)*0.5,] ,pch=21,size=2,aes(x = group,y= gci_count/posi_sites) )+
  facet_wrap(~from,scales = "free_x") +  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=25, family="economica"),
    legend.text =  element_text(color="black", size=25, family="economica"),
    strip.text = element_text(family = "economica", size = 14)
  ) + scale_fill_manual("Gene length (kb)",values = set_color)+ scale_color_manual("Gene length (kb)",values = set_color) + ggtitle(species) +
  scale_x_continuous(breaks=c(-1,-2,-3,-4,-5,1,2,3,4,5),labels = c(0.01,0.1,1,10,100,0.01,0.1,1,10,100)) +
  # scale_x_log10(breaks=c(10,100,1000,10000),labels = c(10,100,"1kb","10kb")) +
  # annotation_logticks(sides="b") +
  annotation_logticks(sides="b") + xlab("Genome distance from start codon or stop codon (log10 scale)") + ylab ("GCi")+ ylim(c(0.3,0.7))
print( p4B )

jpeg(paste(path_pannel,"p4B_supp.jpg",sep=""), width = 5500/1, height = 2500/1,res=400/1)
print(p4B)
dev.off()



############## Supplementary Figure 4

imgA = load.image(paste(path_pannel,"p4A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p4B_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure4_supp.pdf",sep=""), width=5, height=4)
  
  m = matrix(rep(c(1,2),1*2), nrow=2)
  
  layout(m)
  m
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  par(mar=c(0, 2, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B", side=2,at=111, line=1, font=2, cex=1,las=2)
  dev.off()
}
