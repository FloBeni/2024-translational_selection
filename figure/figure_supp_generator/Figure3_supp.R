source("figure/figure_supp generator/library_path.R")


############## Supplementary Pannel 3 A
data12 = read.delim("data/data12.tab")
data12$clade_group = clade_dt[data12$species,]$clade_group

data12 = data12[data12$type_aa == "Wb_WC_notambiguous",]

arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label
data12 = data12[data12$species %in% list_species,]

data12 = data12[ data12$nb_codon_not_decoded == 0 & data12$nb_genes > 5000 & data12$pval_aa_fpkm < 0.05 ,]

data12$ecart = (data12$optifreq_top5 - data12$opti_freq_low50)

p3A = ggplot(data12,aes(y=ecart*100,x=clade_group,fill=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) + 
  geom_boxplot(alpha=.1) + 
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0,angle = 50, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + 
  ylab("Difference in proportion of optimal codons between\nthe top 5% and bottom 50% expressed (%)")  + xlab("") 

p3A = ggMarginal(p3A, type="histogram",fill=set_color[1]) 
p3A
jpeg(paste(path_pannel,"p3A_supp.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p3A)
dev.off()




############## Supplementary Figure 3

imgA = load.image(paste(path_pannel,"p3A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure3_supp.pdf",sep=""), width=6.75, height=4)
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  # mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  dev.off()
}
