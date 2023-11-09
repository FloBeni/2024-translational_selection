source("figure/figure_main generator/library_path.R")


############## Pannel 2 A
data1 = read.delim("data/data1.tab")
dt_graph = data1[data1$species == "Caenorhabditis_elegans", ]

spearman_method_aa = cor.test( dt_graph$tRNASE_copies, dt_graph$obs_codon,method="spearman",exact=F)

p2A = ggplot(dt_graph,aes(x= obs_codon / sum(obs_codon)*100,y=tRNASE_copies,label=amino_acid)) +
  geom_smooth(formula = y ~ x, method="lm", size=1 , col=set_color[1],se=F,linetype='dashed') +
  geom_point(pch=21,size=5,fill=set_color[2]) +
  geom_text(nudge_x = .35,size=6,family="economica") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=20,family="economica"),
    axis.title.y = element_text(color="black", size=20, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
  ) + xlab(paste("Frequency of amino acid weigthed FPKM (%)")) + ylab("tRNA gene copy number") +
  ggtitle(paste(
    "rho = " ,round(spearman_method_aa$estimate,3),
    ", p-value = " ,formatC(spearman_method_aa$p.value, format = "e", digits = 1), sep=""))


jpeg(paste(path_pannel,"p2A.jpg",sep=""), width = 4000/1, height = 2500/1,res=450/1)
print(p2A)
dev.off()


############## Pannel 2 B
dt_graph = data1[!duplicated(data1$species) & data1$species %in% clade_dt$species,]
dt_graph$clade_group = clade_dt[dt_graph$species,]$clade_group

# dt_graph = dt_graph[ dt_graph$nb_gene >= 5000 , ]

p2B = ggplot(dt_graph,aes(y=rho_aa_fpkm,fill=clade_group,x=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = min(dt_graph[dt_graph$rho_aa_fpkm & dt_graph$pval_aa_fpkm < 0.05,]$rho_aa_fpkm) ) +
  geom_text (label="p-value < 0.05", y=.4,x="Lepido Diptera",size=5,family="economica",col="red") + geom_boxplot(alpha=.1) + 
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=1,hjust=1, size=20,angle = 30, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none')+ scale_fill_manual(values=Clade_color) + ylab("Spearmann Rho")  + xlab("") + ylim(0,1)

p2B = ggMarginal(p2B, type="histogram",fill=set_color[1]) 
p2B

jpeg(paste(path_pannel,"p2B.jpg",sep=""), width = 4000/1, height = 2500/1,res=400/1)
print(p2B)
dev.off()

############## Figure 2

imgA = load.image(paste(path_pannel,"p2A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p2B.jpg",sep="") )
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure2.pdf",sep=""), width=4, height=5)
  
  m=matrix(rep(NA,10*10), nrow=10)
  for(i in 1:10){
    m[,i]=c(rep(1,5),rep(2,5))
  }
  layout(m)
  
  par(mar=c(0, 1, 1, 1))
  plot(imgA, axes=FALSE)
  mtext("A", adj=-0.1, side=2,at=1, line=1, font=2, cex=1.2,las=2)
  xcel=400
  ycel=-340
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/1.5-ycel, xright=1000/1.4+xcel, ytop=0-ycel)
  
  par(mar=c(0, 1, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B",adj=-0.1, side=2,at=1, line=1, font=2, cex=1.2,las=2)
  dev.off()
}

