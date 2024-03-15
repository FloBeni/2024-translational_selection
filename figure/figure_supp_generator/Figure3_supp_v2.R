# Generate Supplementary Figure 2
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 2 A

data1_supp = read.delim("data/data7_supp.tab")
dt_graph = data1_supp[data1_supp$species == "Drosophila_melanogaster",]

spearman_method_aa = cor.test( dt_graph$prop_abundance_average, dt_graph$prop_transcriptome_count,method="spearman",exact=F)

pA=ggplot(dt_graph,aes(x = prop_transcriptome_count*100 , y = prop_abundance_average*100)) + geom_smooth(method='lm',linetype='dashed',se=F, formula= y~x,col=set_color[3],size=2) +
  geom_point(pch=21,fill=set_color[4],size=4) + 
  ylab("tRNA abundance (%)") + xlab("Amino-acid frequency (%)") + 
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.55, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  ) + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  )
pA

jpeg(paste(path_pannel,"p2A_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=550/1)
print(pA)
dev.off()


# Supplementary Pannel 2 B

dt_graph = data1_supp[data1_supp$species == "Homo_sapiens",]
spearman_method_aa = cor.test( dt_graph$prop_abundance_average, dt_graph$prop_transcriptome_count,method="spearman",exact=F)

pB=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=prop_abundance_average*100)) + geom_smooth(method='lm',linetype='dashed',se=F, formula= y~x,col=set_color[1],size=2) +
  geom_point(pch=21,fill=set_color[2],size=4) + 
  ylab("tRNA abundance (%)")+ xlab("Amino-acid frequency (%)") + 
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.55, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  ) + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  )
pB

jpeg(paste(path_pannel,"p2B_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=550/1)
print(pB)
dev.off()


# Supplementary Pannel 2 C

dt_graph = data1_supp[data1_supp$species == "Drosophila_melanogaster",]
spearman_method_aa = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)

pC=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=gene_copies)) + geom_smooth(method='lm',linetype='dashed',se=F, formula= y~x,col=set_color[3],size=2) +
  geom_point(pch=21,fill=set_color[4],size=4) + 
  xlab("Amino-acid frequency (%)") + ylab("Gene copy number") + 
   theme_bw() + theme(
     axis.title.x = element_text(color="black", size=26,family="economica"),
     axis.title.y = element_text(color="black", size=26, family="economica"),
     axis.text.y =  element_text(color="black", size=20, family="economica"),
     axis.text.x =  element_text(color="black", size=20, family="economica"),
     title =  element_text(color="black", size=20, family="economica"),
     text =  element_text(color="black", size=31, family="economica"),
     legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
     plot.caption = element_text(hjust = 0.55, face= "italic", size=20, family="economica"),
     plot.caption.position =  "plot"
   ) + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  )
pC


jpeg(paste(path_pannel,"p2C_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=550/1)
print(pC)
dev.off()


# Supplementary Pannel 2 D

dt_graph = data1_supp[data1_supp$species == "Homo_sapiens",]
spearman_method_aa = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)

pD=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=gene_copies)) + geom_smooth(method='lm',linetype='dashed',se=F, formula= y~x,col=set_color[1],size=2) +
  geom_point(pch=21,fill=set_color[2],size=4) + 
  xlab("Amino-acid frequency (%)") + ylab("Gene copy number") + 
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .55, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  ) + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  )
pD

jpeg(paste(path_pannel,"p2D_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=550/1)
print(pD)
dev.off()


# Supplementary Figure 2

imgA = load.image(paste(path_pannel,"p2A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p2B_supp.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p2C_supp.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p2D_supp.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
fly<-readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure3_supp.pdf",sep=""), width=6.75, height=4)
  
  m=matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:10){
      m[i,]=c(rep(1,5),rep(2,5))
    }
  
  for(i in 6:10){
    m[i,]=c(rep(3,5),rep(4,5))
  }
  m
  layout(m)
  
  par(mar=c(0, 1, 1.5, 0))
  plot(imgA, axes=FALSE)
  mtext("A", adj=-.2, side=2,at=-40, line=1, font=2, cex=1.2,las=2)
  xaxis=600/1.1
  yaxis=500/1.1
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1100/4+xaxis, ytop=-900/4+yaxis)
  
  par(mar=c(0, 1, 1.5, 0))
  plot(imgB, axes=FALSE)
  mtext("B", adj=-.2, side=2,at=-40, line=1, font=2, cex=1.2,las=2)
  xhuman=530
  yhuman=-90
  rasterImage(human,xleft=0+xhuman, ybottom=450/1-yhuman, xright=190/1+xhuman, ytop=0-yhuman)

  par(mar=c(0, 1, 1.5, 0))
  plot(imgC, axes=FALSE)
  mtext("C", adj=-.2, side=2,at=-40, line=1, font=2, cex=1.2,las=2)
  xaxis=500/1.1
  yaxis=500/1.1
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1100/4+xaxis, ytop=-900/4+yaxis)

  par(mar=c(0, 1, 1.5, 0))
  plot(imgD, axes=FALSE)
  mtext("D", adj=-.2, side=2,at=-40, line=1, font=2, cex=1.2,las=2)
  xhuman=450
  yhuman=-90
  rasterImage(human,xleft=0+xhuman, ybottom=450/1-yhuman, xright=190/1+xhuman, ytop=0-yhuman)
  dev.off()
}
