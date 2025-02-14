# Generate Supplementary Figure 3
source("figure/figure_supp_generator/library_path.R")
resolution = 2

# Pannel A

data9 = read.delim("data/data9_supp.tab",comment.char = "#")
dt_graph = data9

spearman_method_aa = cor.test( dt_graph$nb_copy_Caenorhabditis_elegans, dt_graph$nb_copy_Hydra_vulgaris,method="spearman",exact=F)

pA = ggplot(dt_graph,aes(x=nb_copy_Hydra_vulgaris ,y=nb_copy_Caenorhabditis_elegans,label=aa_name)) +
  geom_point(pch=21,size=4,fill=set_color[2]) +
  geom_text(nudge_x = 1.7,size=5,family="ubuntu condensed") + theme_bw() +  theme(
    axis.title.x = element_text(color="black", size=22,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=22,hjust=1, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=22, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=22, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=16, family="ubuntu condensed"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.55, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  ) +  
  xlab(expression(paste("tRNA gene copy number (",italic("Hydra vulgaris"),")"))) +
  ylab(expression(paste("tRNA gene copy number (",italic("Caenorhabditis elegans"),")"))) +
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", ", italic("P"), "-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  ) 
pA

jpeg(paste(path_pannel,"p3A_supp.jpg",sep=""), width = 4000/resolution, height = 4000/resolution,res=700/resolution)
print(pA)
dev.off()


# Pannel B

data8 = read.delim("data/data8_supp.tab",comment.char = "#")
data8$clade_group = GTDrift_list_species[data8$species2,]$clade_group

data1 = read.delim("data/data1_supp.tab",comment.char = "#")
rownames(data1) = data1$species
data1 = data1[ data1$nb_genes_filtered >= 5000 ,]
data8$pval_aa_fpkm = data1[data8$species2,]$pval_aa_fpkm < 0.05

dt_graph = data8[data8$species2 %in% data1$species& data8$species2 != "Hydra_vulgaris" ,]

pB = ggplot(dt_graph,aes(y=rho,fill=clade_group,x=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = min(dt_graph[dt_graph$rho & dt_graph$pval < 0.05,]$rho) ) +
  geom_boxplot(alpha=.1) + 
  geom_point(aes(fill=clade_group,pch=pval_aa_fpkm),size=3,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=2, size=22, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=20, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black",vjust=1,hjust=1, size=18,angle = 30, family="ubuntu condensed"),
    title =  element_text(color="black", size=0, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + scale_shape_manual(values=c(24,21)) + ylab("Spearman's rho")  + xlab("") + ylim(0,1)

pB = ggMarginal(pB, type="histogram",fill=set_color[1]) 
pB

jpeg(paste(path_pannel,"p3B_supp.jpg",sep=""), width = 4000/resolution, height = 2500/resolution,res=400/resolution)
print(pB)
dev.off()


# Supplementary Figure 3

imgA = load.image(paste(path_pannel,"p3A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p3B_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure3_supp.pdf",sep=""), width=11, height=4)
  
  m=matrix(rep(NA,10*100), nrow=10)
  for(i in 1:10){
    m[i,]=c(rep(1,35),rep(2,65))
  }
  m
  layout(m)
  
  par(mar=c(0, 2, 2, 1))
  plot(imgA, axes=FALSE)
  mtext("A", adj=0, side=2,at=-200/2, line=1, font=2, cex=2.2,las=2)
  
  par(mar=c(0, 0, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B",adj=-1, side=2,at=150/2, line=1, font=2, cex=2.2,las=2)
  dev.off()
}

