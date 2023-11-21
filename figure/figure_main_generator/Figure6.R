source("figure/figure_main generator/library_path.R")

############## Pannel 6 A
data12 = read.delim("data/data12.tab")
data12$clade_group = clade_dt[data12$species,]$clade_group

table(data12[!duplicated(data12$species),]$nb_codon_not_decoded)
data12 = data12[ data12$nb_codon_not_decoded == 0 & data12$nb_genes > 5000 & data12$pval_aa_fpkm < 0.05 ,]

data12$ecart = (data12$optifreq_top5-data12$opti_freq_low50) - (data12$optifreq_top5_intron-data12$opti_freq_low50_intron)

data4 = read.delim("data/data4.tab")
data4$clade_group = clade_dt[data4$species,]$clade_group
data4 = data4[data4$species %in% data12$species ,]


dt_graph = data4[ grepl("WC_duet_ambiguous",data4$type_aa) ,]
dt_graph = dt_graph[ !duplicated(dt_graph$species) ,]

dt_graph$prop_GU = dt_graph$prop_GU  / dt_graph$nb_aa

p6A = ggplot(dt_graph,aes(x=prop_GU))  + ggtitle("WC_duet_ambiguous") +
  geom_histogram(col="black",fill=set_color[1]) + 
  # facet_wrap(~clade_group) + 
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(size = 20, family = "economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color)+ylab("Proportion of ambiguous duet with wobble GU")  + xlab("") +
  geom_text(x=0.68,y=15,label="Asterias rubens",size = 10, family = "economica")

p6A

jpeg(paste(path_pannel,"p6A.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p6A)
dev.off()


############## Pannel 6 B
dt_graph = data12[data12$type_aa == "WC_duet_ambiguous",]


p6B = ggplot(dt_graph,aes(y=ecart*100,x=clade_group,fill=clade_group))  + ggtitle("WC_duet_ambiguous") +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) + 
  geom_boxplot(alpha=.1) + 
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + 
  ylab("Difference in proportion of WCP codons between\nthe top 5% and bottom 50% expressed")  + xlab("") 

p6B = ggMarginal(p6B, type="histogram",fill=set_color[1]) 
p6B


jpeg(paste(path_pannel,"p6B.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p6B)
dev.off()


############## Pannel 6 C

data6 = read.delim("data/data6.tab")
data6$clade_group = clade_dt[data6$species,]$clade_group

data6 = data6[data6$species %in% data12$species,]
data6 = data6[data6$nb_gene > 500,]

data6$ecart = data6$mean_highconstrain - data6$mean_unconstrain

dt_graph = data6[data6$type_aa == "WC_duet_ambiguous",]


p6C = ggplot(dt_graph,aes(y=ecart*100,x=clade_group,fill=clade_group))  + ggtitle("WC_duet_ambiguous") +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color)+  
  ylab("Difference in proportion of WCP codons between\nthe 25% highest and the 25% lowest constrained sites") + xlab("") 

p6C = ggMarginal(p6C, type="histogram",fill=set_color[1])
p6C

jpeg(paste(path_pannel,"p6C.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p6C)
dev.off()

############## Pannel 6 D
dt_graph = data4[ grepl("WB_vs_WC_notduet",data4$type_aa) ,]
dt_graph = dt_graph[ !duplicated(dt_graph$species) ,]

dt_graph$prop_GU = dt_graph$prop_GU  / dt_graph$nb_aa

p6D = ggplot(dt_graph,aes(x=prop_GU,fill=clade_group))  + ggtitle("WB_vs_WC_notduet") +
  geom_histogram(col="black") + facet_wrap(~clade_group) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=10, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica"),
    strip.text = element_text(size = 20, family = "economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color)+ylab("Proportion of ambiguous (diff duet with wobble IC")  + xlab("") 

p6D

jpeg(paste(path_pannel,"p6D.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p6D)
dev.off()


############## Pannel 6 E
dt_graph = data12[data12$type_aa == "IC",]

p6E = ggplot(dt_graph,aes(y=ecart*100,x=clade_group,fill=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) + 
  geom_boxplot(alpha=.1) + ggtitle("IC") +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + 
  ylab("Difference in proportion of optimal codons between\nthe top 5% and bottom 50% expressed") + xlab("") 

p6E = ggMarginal(p6E, type="histogram",fill=set_color[1]) 
p6E


jpeg(paste(path_pannel,"p6E.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p6E)
dev.off()


############## Pannel 6 F
dt_graph = data6[data6$type_aa == "IC",]

p6F = ggplot(dt_graph,aes(y=ecart*100,x=clade_group,fill=clade_group))  + ggtitle("IC") +
  geom_hline(size=1,linetype="dashed",col="red",
             yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color)+
  ylab("Difference in proportion of optimal codons between\nthe 25% highest and the 25% lowest constrained sites")+ xlab("")

p6F = ggMarginal(p6F, type="histogram",fill=set_color[1])
p6F

jpeg(paste(path_pannel,"p6F.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(p6F)
dev.off()



############## Figure 6

imgA = load.image(paste(path_pannel,"p6A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p6B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p6C.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p6D.jpg",sep="") )
imgE = load.image(paste(path_pannel,"p6E.jpg",sep="") )
imgF = load.image(paste(path_pannel,"p6F.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure6.pdf",sep=""), width=10, height=7)
  
  m = matrix(rep(NA,20*15), nrow=15)
  
  for(i in 1:10){
    m[,i]=c(rep(1,5),rep(2,5),rep(3,5))
  }
  for(i in 11:20){
    m[,i]=c(rep(4,5),rep(5,5),rep(6,5))
  }
  layout(m)
  m
  
  par(mar=c(1, 0, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=100,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgD, axes=FALSE)
  mtext("D",at=100,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgE, axes=FALSE)
  mtext("E",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 0, 0))
  plot(imgF, axes=FALSE)
  mtext("F",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  dev.off()
}

