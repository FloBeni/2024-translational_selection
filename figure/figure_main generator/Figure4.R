source("figure/figure_main generator/library_path.R")

############## Pannel 3 A
data4 = read.delim("data/data4.tab")
data4$set = str_replace_all(data4$set,"all,","genes,")
dt_graph = data4[ data4$species == "Homo_sapiens" & grepl("Wb_WC_notambiguous",data4$type_aa) ,]
dt_graph = data4[ data4$species == "Musca_domestica" & grepl("Wb_WC_notambiguous",data4$type_aa) ,]


p4A = ggplot(data4 ,
             aes(x=fpkm ,y=freq*100,linetype=set,fill=categorie,col=categorie))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2) +
  geom_point(data=dt_graph,aes(pch=set),col="black",size=4)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Frequency optimal codons (%)") + theme_bw() + theme(
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
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("optimal codons",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#E31A1C") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("intronic",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#FB9A99") +
  geom_point(data =  dt_graph[grepl("optimal codons",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#E31A1C",size=6)+
  geom_point(data = dt_graph[grepl("intronic",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#FB9A99",size=6)+
  ggtitle(unique(dt_graph$set)) +   guides(fill="none",color="none",linetype="none",shape="none")+ annotation_logticks(sides = "b")

p4A

jpeg(paste(path_pannel,"p4A.jpg",sep=""),  width = 8000/2,  5400/2,res=900/2)
print(p4A)
dev.off()

############## Pannel 3 B

dt_graph = data4[ data4$species == "Caenorhabditis_elegans" & grepl("Wb_WC_notambiguous",data4$type_aa) ,]
# dt_graph = data4[ data4$species == "Drosophila_melanogaster" & grepl("Wb_WC_notambiguous",data4$type_aa) ,]

p4B = ggplot(data4 ,
             aes(x=fpkm ,y=freq*100,linetype=set,fill=categorie,col=categorie))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2) +
  geom_point(data=dt_graph,aes(pch=set),col="black",size=4)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Frequency optimal codons (%)") + theme_bw() + theme(
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
  )+          scale_x_log10(
    breaks=c(0.005,0.01,0.1,1,10,100,500,1000,10000,50000),
    labels=c(0.005,0.01,0.1,1,10,100,500,1000,10000,50000))+ ylim(0.3*100,0.7*100) +
  # geom_hline(yintercept = 0.48283,size=1,linetype="dashed",col="#E31A1C") +
  # geom_hline(yintercept = 0.40359,size=1,linetype="dashed",col="#FB9A99") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("optimal codons",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#E31A1C") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("intronic",dt_graph$categorie),]$freq)*100,size=1,linetype="dashed",col="#FB9A99") +
  geom_point(data =  dt_graph[grepl("optimal codons",dt_graph$categorie) &dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#E31A1C",size=6)+
  geom_point(data = dt_graph[grepl("intronic",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) , ],col="black",pch=21,fill="#FB9A99",size=6) +
  ggtitle(unique(dt_graph$set)) +   guides(linetype="none",shape="none")+ annotation_logticks(sides = "b")

jpeg(paste(path_pannel,"p4B.jpg",sep=""),  width = 11000/2,  5500/2,res=900/2)
print(p4B)
dev.off()

############## Pannel 3 C

data12 = read.delim("data/data12.tab")
data12$clade_group = clade_dt[data12$species,]$clade_group

data12 = data12[data12$type_aa == "Wb_WC_notambiguous",]

data12 = data12[ data12$nb_codon_not_decoded == 0 & data12$pval_aa_fpkm < 0.05 ,]

data12$ecart = (data12$optifreq_top5-data12$opti_freq_low50) - (data12$optifreq_top5_intron-data12$opti_freq_low50_intron)

p4C = ggplot(data12,aes(y=ecart*100,x=clade_group,fill=clade_group,label=species))  +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) + 
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0,angle = 50, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + 
  ylab("Difference in proportion of optimal codons between\nthe top 5% and bottom 50% expressed (%)")  + xlab("")  + ylim(-5,20)

p4C = ggMarginal(p4C, type="histogram",fill=set_color[1]) 
p4C

jpeg(paste(path_pannel,"p4C.jpg",sep=""), width = 5500/1, height = 3000/1,res=500/1)
print(p4C)
dev.off()




############## Figure 3

imgA = load.image(paste(path_pannel,"p4A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p4B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p4C.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))
Drosophila_melanogaster = readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))
clade_png<-readPNG(paste(path_require,"clade.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure4.pdf",sep=""), width=10, height=6)
  
  m = matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:5){
    m[i,]=c(rep(1,5),rep(2,5))
  }
  
  for(i in 6:10){
    m[i,]=c(rep(3,10))
  }
  layout(m)
  
  par(mar=c(0, 5, 0, 5))
  plot(imgA, axes=FALSE)
  mtext("A",at=100,adj=1, side=2, line=1, font=2, cex=1.7,las=2)
  xhuman=450
  yhuman=-300
  rasterImage(human,xleft=0+xhuman, ybottom=350/1.1-yhuman, xright=190/1.1+xhuman, ytop=0-yhuman)
  par(mar=c(0, 0, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=100,adj=1, side=2, line=1, font=2, cex=1.7,las=2)
  xcel=400
  ycel=-350
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/2-ycel, xright=1000/2+xcel, ytop=0-ycel)
  # xcel=100
  # ycel=-140
  # rasterImage(Drosophila_melanogaster,xleft=0+xcel, ybottom=350/3.5-ycel, xright=450/3.5+xcel, ytop=0-ycel)
  par(mar=c(0, 0, 0, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=100,adj=-10, side=2, line=1, font=2, cex=1.7,las=2)
  xmonkey=5000
  ymonkey=300
  rasterImage(clade_png,xleft=0+xmonkey, ybottom=800/0.5+ymonkey, xright=400/0.5+xmonkey, ytop=ymonkey)
  dev.off()
}

