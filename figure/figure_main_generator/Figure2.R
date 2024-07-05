# Generate Figure 2
source("figure/figure_main_generator/library_path.R")


# Pannel A

data3 = read.delim("data/data3_supp.tab",comment.char = "#")
dt_graph = data3[data3$species == "Caenorhabditis_elegans", ]

spearman_method_aa = cor.test( dt_graph$tRNA_gene_copy, dt_graph$obs_aminoacid,method="spearman",exact=F)

pA = ggplot(dt_graph,aes(x= obs_aminoacid / sum(obs_aminoacid) * 100,y=tRNA_gene_copy,label=amino_acid)) +
  geom_point(pch=21,size=4,fill=set_color[2]) +
  geom_text(nudge_x = .35,size=5,family="ubuntu condensed") + theme_bw() +  theme(
    axis.title.x = element_text(color="black",vjust=0, size=26,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=25, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=22, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=22, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=16, family="ubuntu condensed"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.55,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  ) + xlab(paste("Amino-acid frequency (%)")) + ylab("tRNA gene copy number") + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  ) 
pA

jpeg(paste(path_pannel,"p2A.jpg",sep=""), width = 4000/2, height = 4000/2,res=700/2)
print(pA)
dev.off()


# Pannel B

data1 = read.delim("data/data1_supp.tab",comment.char = "#")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_genes_filtered >= 5000,]

dt_graph = data1

pB = ggplot(dt_graph,aes(y=rho_aa_fpkm,fill=clade_group,x=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = min(dt_graph[dt_graph$rho_aa_fpkm & dt_graph$pval_aa_fpkm < 0.05,]$rho_aa_fpkm) ) +
  geom_text (label="p-value < 0.05", y=.4,x="Lepido Diptera",size=5,family="ubuntu condensed",col="red") + geom_boxplot(alpha=.1) + 
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=22, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=20, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black",vjust=1,hjust=1, size=18,angle = 30, family="ubuntu condensed"),
    title =  element_text(color="black", size=0, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + ylab("Spearmann rho")  + xlab("") + ylim(0,1)

pB = ggMarginal(pB, type="histogram",fill=set_color[1]) 
pB

jpeg(paste(path_pannel,"p2B.jpg",sep=""), width = 4000/2, height = 2500/2,res=400/2)
print(pB)
dev.off()


# Figure 2

imgA = load.image(paste(path_pannel,"p2A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p2B.jpg",sep="") )
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure2.pdf",sep=""), width=11, height=4)
  
  m=matrix(rep(NA,10*100), nrow=10)
  for(i in 1:10){
    m[i,]=c(rep(1,35),rep(2,65))
  }
  m
  layout(m)
  
  par(mar=c(0, 2, 2, 1))
  plot(imgA, axes=FALSE)
  mtext("A", adj=0, side=2,at=-200/2, line=1, font=2, cex=2.2,las=2)
  xcel=630/2
  ycel=-140
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/1.5/2-ycel, xright=1000/1.4/2+xcel, ytop=0-ycel)
  
  par(mar=c(0, 0, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B",adj=-1, side=2,at=150/2, line=1, font=2, cex=2.2,las=2)
  dev.off()
}

