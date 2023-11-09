source("figure/figure_supp generator/library_path.R")


############## Supplementary Pannel 5 A
arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label

# Freq_opti = read.delim("data/Freq_opti.tab")
# Freq_opti = Freq_opti[Freq_opti$species %in% list_species,]
# Freq_opti$clade_group = clade_dt[Freq_opti$species,]$clade_group
# Freq_opti = Freq_opti[ Freq_opti$nb_codon_not_decoded == 0 & Freq_opti$nb_genes > 5000 & Freq_opti$pval_aa_fpkm < 0.05 ,]

data10 = read.delim(paste("data/data10.tab",sep="") , header=T )
data10 = read.delim(paste("/home/fbenitiere/data/Projet-NeGA/translational_selection/GC_gap_per_window_per_species_100bp_bins_all_genes.tab",sep="") , header=T )

data10$clade_group = clade_dt[data10$species,]$clade_group
data10 = data10[data10$species %in% list_species ,]
unique(data10$species)
data10$from = str_replace_all(data10$from , "from_5prime","From start")
data10$from = str_replace_all(data10$from , "from_3prime","From stop")
data10$from = factor(data10$from,levels = c("From start" , "From stop"))

data10$label_facet = paste(data10$clade_group, ", N=",  table(data10[!duplicated(data10$species),]$clade_group)[data10$clade_group])
vector_label = unique(data10$label_facet)
names(vector_label) = unique(data10$clade_group)

data10$label_facet = factor(data10$label_facet , levels = vector_label[levels(data10$clade_group)] )

# dg = data10[data10$clade_group == "Mammalia" & data10$nb_genes >= 200, ]
dg = data10[ data10$nb_genes >= 200, ]
dg = data10[ data10$nb_genes >= 1000, ]
dg$ratio  = dg$gci_count/dg$posi_sites

data10[data10$from == "From start",]$group = log10(data10[data10$from == "From start",]$group)
data10[data10$from == "From stop",]$group = -log10(data10[data10$from == "From stop",]$group)

p5A = ggplot(data10[data10$nb_genes >= 1000,],aes(x = group,group=clade_group,y= gci_count/posi_sites  , label=nb_genes,fill=clade_group,col=clade_group)) +
  geom_line(lwd = 0.5,aes(group=species), alpha=0.7) + 
  # geom_line(data = data10[data10$nb_genes >= 200 & data10$species == "Thrips_palmi",] , lwd = 0.5,aes(group=species), alpha=0.7,col="red") + 
  
  
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
  # scale_x_log10(breaks=c(10,100,1000,10000),labels = c(10,0.1,"1","10")) +
  # annotation_logticks(sides="b") +
  xlab("Genomic distance from start codon or stop codon (kb, log scale)") + ylab ("GCi")+ theme(legend.position='none') 
print(p5A)


jpeg(paste(path_pannel,"p5A_supp.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p5A)
dev.off()



############## Supplementary Figure 5

imgA = load.image(paste(path_pannel,"p5A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure5_supp.pdf",sep=""), width=6.75, height=4)
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  # mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  dev.off()
}