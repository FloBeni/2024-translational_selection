source("figure/figure_supp generator/library_path.R")


############## Supplementary Pannel 2 A
data_supp_1 = read.delim("data/data_supp_1.tab")
dt_graph = data_supp_1[data_supp_1$species == "Drosophila_melanogaster",]

p2A=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=prop_abundance_average*100)) + geom_smooth(method='lm', formula= y~x,col=set_color[3],size=2) +
  geom_point(pch=21,fill=set_color[4],size=3) + 
  ylab("Proportion uniquely aligned reads (%)")+ xlab("Frequency of amino acid weigthed FPKM (%)") + 
  geom_abline() + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=21,family="serif"),
    axis.title.y = element_text(color="black", size=21, family="serif"),
    axis.text.y =  element_text(color="black", size=20, family="serif"),
    axis.text.x =  element_text(color="black", size=20,hjust=1, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.text =  element_text(color="black", size=18, family="serif")
  ) +   ggtitle(paste("D .mel, Rho Spearman =",round(cor.test( dt_graph$prop_transcriptome_count, dt_graph$prop_abundance_average,method="spearman",exact=F)$estimate,2)))
p2A

jpeg(paste(path_pannel,"p2A_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=500/1)
print(p2A)
dev.off()


############## Supplementary Pannel 2 B
dt_graph = data_supp_1[data_supp_1$species == "Homo_sapiens",]

p2B=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=prop_abundance_average*100)) + geom_smooth(method='lm', formula= y~x,col=set_color[1],size=2) +
  geom_point(pch=21,fill=set_color[2],size=3) + 
  ylab("Proportion uniquely aligned reads (%)")+ xlab("Frequency of amino acid weigthed FPKM (%)") + 
  geom_abline() + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=21,family="serif"),
    axis.title.y = element_text(color="black", size=21, family="serif"),
    axis.text.y =  element_text(color="black", size=20, family="serif"),
    axis.text.x =  element_text(color="black", size=20,hjust=1, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.text =  element_text(color="black", size=18, family="serif")
  ) +   ggtitle(paste("Homo .sap, Rho Spearman =",round(cor.test( dt_graph$prop_transcriptome_count, dt_graph$prop_abundance_average,method="spearman",exact=F)$estimate,2)))
p2B

jpeg(paste(path_pannel,"p2B_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=500/1)
print(p2B)
dev.off()

############## Supplementary Pannel 2 C
dt_graph = data_supp_1[data_supp_1$species == "Drosophila_melanogaster",]

p2C=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=gene_copies)) + geom_smooth(method='lm', formula= y~x,col=set_color[3],size=2) +
  geom_point(pch=21,fill=set_color[4],size=3) + 
  xlab("Frequency of amino acid weigthed FPKM (%)") + ylab("Gene copy number") + 
   theme_bw() + theme(
    axis.title.x = element_text(color="black", size=21,family="serif"),
    axis.title.y = element_text(color="black", size=21, family="serif"),
    axis.text.y =  element_text(color="black", size=20, family="serif"),
    axis.text.x =  element_text(color="black", size=20,hjust=1, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.text =  element_text(color="black", size=18, family="serif")
  ) +   ggtitle(paste("D .mel, Rho Spearman =",round(cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)$estimate,2))) +
  ylim(0,50)
p2C


jpeg(paste(path_pannel,"p2C_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=500/1)
print(p2C)
dev.off()

############## Supplementary Pannel 2 D
dt_graph = data_supp_1[data_supp_1$species == "Homo_sapiens",]

p2D=ggplot(dt_graph,aes(x=prop_transcriptome_count*100,y=gene_copies)) + geom_smooth(method='lm', formula= y~x,col=set_color[1],size=2) +
  geom_point(pch=21,fill=set_color[2],size=3) + 
  xlab("Frequency of amino acid weigthed FPKM (%)") + ylab("Gene copy number") + 
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=21,family="serif"),
    axis.title.y = element_text(color="black", size=21, family="serif"),
    axis.text.y =  element_text(color="black", size=20, family="serif"),
    axis.text.x =  element_text(color="black", size=20,hjust=1, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.text =  element_text(color="black", size=18, family="serif")
  ) +   ggtitle(paste("Homo .sap, Rho Spearman =",round(cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)$estimate,2))) +
  ylim(0,50)
p2D

jpeg(paste(path_pannel,"p2D_supp.jpg",sep=""), width = 4000/1, height = 2500/1,res=500/1)
print(p2D)
dev.off()


############## Supplementary Figure 2

imgA = load.image(paste(path_pannel,"p2A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p2B_supp.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p2C_supp.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p2D_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure2_supp.pdf",sep=""), width=6.75, height=4)
  
  m=matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:10){
      m[i,]=c(rep(1,5),rep(2,5))
    }
  
  for(i in 6:10){
    m[i,]=c(rep(3,5),rep(4,5))
  }
  m
  layout(m)
  
  par(mar=c(0, 2, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  
  par(mar=c(0, 2, 2, 0))
  plot(imgB, axes=FALSE)
  mtext("B", side=2,at=111, line=1, font=2, cex=1,las=2)
  
  par(mar=c(0, 2, 2, 0))
  plot(imgC, axes=FALSE)
  mtext("C", side=2,at=1, line=1, font=2, cex=1,las=2)
  
  par(mar=c(0, 2, 2, 0))
  plot(imgD, axes=FALSE)
  mtext("D", side=2,at=1, line=1, font=2, cex=1,las=2)
  dev.off()
}
