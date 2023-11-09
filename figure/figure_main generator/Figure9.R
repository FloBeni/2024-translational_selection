source("figure/figure_main generator/library_path.R")


############## Pannel 8 A
Freq_opti = read.delim("data/Freq_opti.tab")
Freq_opti[,c("clade_group","lifespan","length","weight")] = clade_dt[Freq_opti$species,c("clade_group","lifespan","length","weight")]

arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label
Freq_opti = Freq_opti[Freq_opti$species %in% list_species,]

dt_graph = Freq_opti[Freq_opti$type_aa == "Wb_WC_notambiguous",]
dt_graph = dt_graph[dt_graph$species %in% clade_dt$species,]
dt_graph = dt_graph[ dt_graph$nb_codon_not_decoded == 0 & dt_graph$nb_genes > 5000 & dt_graph$pval_aa_fpkm < 0.05 ,]
colnames(data7)

code = code[code$nb_syn != 2 | !duplicated(code$aa_name),]
data7 = data7[data7$codon %in% code$codon & data7$sum != 1 , ]
for (codon in unique(data7$codon)){
  p8A = ggplot(data7,aes(y=per_gene,x=per_genei,fill=clade_group))  +
    geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
      axis.title.x = element_text(color="black", size=25,family="economica"),
      axis.title.y = element_text(color="black", size=25, family="economica"),
      axis.text.y =  element_text(color="black", size=20, family="economica"),
      axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
      title =  element_text(color="black", size=15, family="economica"),
      legend.text =  element_text(color="black", size=20, family="economica")
    ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + facet_wrap(~codon,scales = "free")
  print(p8A)
}

(paste(path_pannel,"p8A.jpg",sep=""), width = 3000/1, height = 2500/1,res=400/1)
print(p8A)
dev.off()


############## Pannel 8 B
p8B = ggplot(dt_graph,aes(y=gc3,x=gci,fill=clade_group))  +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + geom_abline()
p8B

jpeg(paste(path_pannel,"p8B.jpg",sep=""), width = 3000/1, height = 2500/1,res=400/1)
print(p8B)
dev.off()


############## Pannel 8 C
data7 = read.delim("data/data14.tab")
data7$clade_group = clade_dt[data7$species,]$clade_group

arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label
data7 = data7[data7$species %in% list_species,]

p8C = ggplot(data7[ data7$genome_character == "CG",],aes(y=shift*100,x=clade_group,fill=clade_group)) +
  geom_hline(yintercept=0,fill="red",alpha=0.2)+ geom_boxplot() + 
  scale_fill_manual("Clades",values=Clade_color) + scale_size(range=c(0,8),guide="legend") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20,angle=60, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + ylab("Shift observed - expected (%)") + xlab("")


jpeg(paste(path_pannel,"p8C.jpg",sep=""), width = 4500/1, height = 2500/1,res=250/1)
print(p8C)
dev.off()


############## Figure 8

imgA = load.image(paste(path_pannel,"p8A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p8B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p8C.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure8.pdf",sep=""), width=4, height=7)
  
  m = matrix(rep(NA,10*15), nrow=15)
  
  for(i in 1:10){
    m[,i]=c(rep(1,5),rep(2,5),rep(3,5))
  }
  layout(m)
  m
  
  par(mar=c(1, 0, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  dev.off()
}

