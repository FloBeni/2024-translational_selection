# Generate Figure 4
source("figure/figure_main_generator/library_path.R")


# Pannel 4 A

data5 = read.delim("data/data5.tab")
data5$set = str_replace_all(data5$set,"all,","genes,")
dt_graph = data5[ data5$species == "Homo_sapiens" & grepl("Wb_WC_notambiguous",data5$type_aa) ,]


pA = ggplot(data5 ,
             aes(x=fpkm ,y=freq*100,fill=categorie,col=categorie))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2) +
  geom_point(data=dt_graph,pch=21,col="black",size=4)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Codon set frequency (%)") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=23, family="economica"),
    axis.text.x =  element_text(color="black", size=23, family="economica"),
    title =  element_text(color="black", size=27, family="economica"),
    legend.text =  element_text(color="black", size=25, family="economica"),
    strip.text = element_text(size=15)
  )+labs(fill='Categories',color='Categories',shape='',linetype='')+
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  ) + scale_x_log10(
    breaks=c(0.005,0.01,0.1,1,10,100,500,1000,10000,50000),
    labels=c(0.005,0.01,0.1,1,10,100,500,1000,10000,50000)) + ylim(0.3*100,0.7*100) +
  # geom_hline(yintercept = 0.50129,size=1,linetype="dashed",col="#E31A1C") +
  # geom_hline(yintercept = 0.47242,size=1,linetype="dashed",col="#FB9A99") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("(POC)",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#E31A1C") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("(POCMT)",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#FB9A99") +
  geom_point(data =  dt_graph[grepl("(POC)",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#E31A1C",size=6)+
  geom_point(data = dt_graph[grepl("(POCMT)",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#FB9A99",size=6)+
  ggtitle(paste(unique(dt_graph$set)," genes",sep="")) +   guides(fill="none",color="none",linetype="none",shape="none")+ annotation_logticks(sides = "b")

pA

jpeg(paste(path_pannel,"p4A.jpg",sep=""),  width = 7500/2,  5400/2,res=900/1.8)
print(pA)
dev.off()


# Pannel 4 B

dt_graph = data5[ data5$species == "Caenorhabditis_elegans" & grepl("Wb_WC_notambiguous",data5$type_aa) ,]

pB = ggplot(data5 ,
             aes(x=fpkm ,y=freq*100,fill=categorie,col=categorie))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2) +
  geom_point(data=dt_graph,pch=21,col="black",size=4)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Codon set frequency (%)") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=23, family="economica"),
    axis.text.x =  element_text(color="black", size=23, family="economica"),
    title =  element_text(color="black", size=27, family="economica"),
    legend.text =  element_text(color="black", size=25, family="economica",vjust = 1.5,margin = margin(t = 10)),
    strip.text = element_text(size=15)
  )+labs(fill='Categories',color='Categories',shape='',linetype='')+
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  )+          scale_x_log10(
    breaks=c(0.005,0.01,0.1,1,10,100,500,1000,10000,50000),
    labels=c(0.005,0.01,0.1,1,10,100,500,1000,10000,50000))+ ylim(0.3*100,0.7*100) +
  # geom_hline(yintercept = 0.48283,size=1,linetype="dashed",col="#E31A1C") +
  # geom_hline(yintercept = 0.40359,size=1,linetype="dashed",col="#FB9A99") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("(POC)",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#E31A1C") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("(POCMT)",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#FB9A99") +
  geom_point(data =  dt_graph[grepl("(POC)",dt_graph$categorie) &dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#E31A1C",size=6)+
  geom_point(data = dt_graph[grepl("(POCMT)",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#FB9A99",size=6) +
  ggtitle(paste(unique(dt_graph$set)," genes",sep="")) +   guides(linetype="none",shape="none")+ annotation_logticks(sides = "b")
pB

jpeg(paste(path_pannel,"p4B.jpg",sep=""),  width = 12000/2,  5500/2,res=900/1.8)
print(pB)
dev.off()


# Pannel 4 C

data1 = read.delim("data/data1.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0 & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000,]

pC = ggplot(data1,aes(y=expressed_overused_background_WB_WC_notambiguous,x=clade_group,fill=clade_group,label=species))  +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica",margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0,angle = 50, family="economica"),
    title =  element_text(color="black", size=0, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Difference in proportion of POC between\nthe top 5% and bottom 50% expressed (%)")  + xlab("")  + ylim(-5,20)

pC = ggMarginal(pC, type="histogram",fill=set_color[1])
pC

jpeg(paste(path_pannel,"p4C.jpg",sep=""), width = 5500/1, height = 3000/1,res=460/1)
print(pC)
dev.off()


# Figure 4

imgA = load.image(paste(path_pannel,"p4A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p4B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p4C.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))
Drosophila_melanogaster = readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))
clade_png<-readPNG(paste(path_require,"clade.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure4.pdf",sep=""), width=9, height=6)
  
  m = matrix(rep(NA,100*10), nrow=10)
  
  for(i in 1:5){
    m[i,]=c(rep(1,45),rep(2,55))
  }
  
  for(i in 5:10){
    m[i,]=c(rep(3,10))
  }
  layout(m)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=100,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)
  xhuman=470
  yhuman=-350
  rasterImage(human,xleft=0+xhuman, ybottom=450/.9-yhuman, xright=190/.9+xhuman, ytop=0-yhuman)
  par(mar=c(0, 0, 1, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=100,adj=0.5, side=2, line=1, font=2, cex=1.7,las=2)
  xcel=500
  ycel=-350
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/1.5-ycel, xright=1000/1.5+xcel, ytop=0-ycel)
  par(mar=c(1,12, 0, 13))
  plot(imgC, axes=FALSE)
  mtext("C",at=200,adj=1, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(1,0, 0, 0))
  xmonkey=5500
  ymonkey=500
  rasterImage(clade_png,xleft=0+xmonkey, ybottom=800/.4+ymonkey, xright=500/.4+xmonkey, ytop=ymonkey)
  dev.off()
}

