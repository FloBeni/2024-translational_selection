# Generate Figure 4
source("figure/figure_main_generator/library_path.R")

resolution=3
# Pannel A

data5 = read.delim("data/data5_supp.tab",comment.char = "#")
data5$gene_set = str_replace_all(data5$gene_set,"all,","genes,")
data5[data5$categorie == "POC-matching triplets (POCMT)",]$categorie = "control"
data5[data5$categorie == "Putative optimal codons (POC)",]$categorie = ""

dt_graph = data5[ data5$species == "Homo_sapiens" & data5$set != "POCs",]

pA = ggplot(dt_graph ,
            aes(x=fpkm ,y=freq*100,fill=paste(set,categorie),col=paste(set,categorie)))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2,aes(linetype=paste(set,categorie))) +
  geom_point(data=dt_graph,pch=21,col="black",size=3)+
  scale_fill_manual(values=set_color[c(2,1,4,3)]) +
  scale_color_manual(values=set_color[c(2,1,4,3)]) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Codon frequency (%)") + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=0, size=20,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=2, size=20, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=20, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=20, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=23, family="ubuntu condensed"),
    strip.text = element_text(size=15)
  )+labs(fill='Categories',color='Categories',shape='',linetype='')+
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  ) + 
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & !grepl("control",dt_graph$categorie) & grepl("POC1",dt_graph$set),]$freq)*100,size=1,linetype="dashed",col="#E31A1C") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("control",dt_graph$categorie) & grepl("POC1",dt_graph$set),]$freq)*100,size=1,linetype="dashed",col="#FB9A99") +
  geom_point(data =  dt_graph[!grepl("control",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) & grepl("POC1",dt_graph$set) , ],col="black",pch=21,fill="#E31A1C",size=6)+
  geom_point(data = dt_graph[grepl("control",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) & grepl("POC1",dt_graph$set), ],col="black",pch=21,fill="#FB9A99",size=6)+
  scale_x_log10(
    breaks=c(0.01,0.1,1,10,100,1000,10000,50000),
    labels=c(0.01,0.1,1,10,100,1000,10000,50000),limits=c(0.005,1000))+ ylim(0.2*100,0.8*100) +
  ggtitle(paste(unique(dt_graph$gene_set)," genes",sep="")) + annotation_logticks(sides = "b")+   guides(fill="none",color="none",linetype="none",shape="none")
pA

jpeg(paste(path_pannel,"p4A.jpg",sep=""),  width = 8500/2/resolution,  5400/2/resolution,res=1000/1.8/resolution)
print(pA)
dev.off()


# Pannel B

dt_graph = data5[ data5$species == "Caenorhabditis_elegans" & data5$set != "POCs",]

data1 = read.delim("data/data1_supp.tab",comment.char = "#")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]

pB = ggplot(dt_graph ,
            aes(x=fpkm ,y=100*freq,fill=paste(set,categorie),col=paste(set,categorie)))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2,aes(linetype=paste(set,categorie))) +
  geom_point(data=dt_graph,pch=21,col="black",size=3)+
  scale_fill_manual(values=set_color[c(2,1,4,3,8)]) +
  scale_color_manual(values=set_color[c(2,1,4,3,8)]) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Codon frequency (%)") + theme_bw() +
  theme(
    axis.title.x = element_text(color="black",vjust=0, size=20,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=2, size=20, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=20, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=20, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    legend.title =  element_text(color="black", size=23, family="ubuntu condensed"),
    strip.text = element_text(size=15)
  )+labs(fill='Categories',color='Categories',shape='',linetype='')+
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  )  + scale_x_log10(
    breaks=c(0.01,0.1,1,10,100,1000,10000,50000),
    labels=c(0.01,0.1,1,10,100,1000,10000,50000),limits=c(0.005,1000))+ ylim(0.2*100,0.8*100) +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & !grepl("control",dt_graph$categorie) & grepl("POC1",dt_graph$set),]$freq)*100,size=1,linetype="dashed",col="#E31A1C") +
  geom_hline(yintercept = mean(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) & grepl("control",dt_graph$categorie) & grepl("POC1",dt_graph$set),]$freq)*100,size=1,linetype="dashed",col="#FB9A99") +
  geom_point(data =  dt_graph[!grepl("control",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) & grepl("POC1",dt_graph$set) , ],col="black",pch=21,fill="#E31A1C",size=6)+
  geom_point(data = dt_graph[grepl("control",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) & grepl("POC1",dt_graph$set), ],col="black",pch=21,fill="#FB9A99",size=6)+
  ggtitle(paste(unique(dt_graph$gene_set)," genes",sep="")) +   guides(linetype="none",shape="none",fill="none")+ annotation_logticks(sides = "b") + 
  guides(color = guide_legend(byrow = TRUE,override.aes = list(lwd=5))) +  theme(legend.spacing.y = unit(.1, 'cm')) 
# geom_line(data = data.frame(Fop_estimate = Fop_estimation(dt_graph$fpkm),fpkm=dt_graph$fpkm),
#           aes(x=fpkm ,y=100*Fop_estimate,fill="simulation",col="simulation"),size=2) 
pB

jpeg(paste(path_pannel,"p4B.jpg",sep=""),  width = 11000/2/resolution,  5500/2/resolution,res=1000/1.8/resolution)
print(pB)
dev.off()


# Pannel C

dt_graph = data1
ylabel = "expressed_overused_background_POC2"
xlabel = "expressed_overused_background_POC1"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,]

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T,lm_obliged=F)

pC =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species"))  +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_abline(linetype="dashed") +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=-0, size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=25, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=26, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=26, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.35, face= "italic", size=23, family="ubuntu condensed"),
    plot.caption.position =  "plot",
    legend.title =  element_text(color="black", size=25, family="ubuntu condensed"),
  )+ guides(fill = guide_legend(override.aes = list(size=5))) +
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", ", italic("P"), "-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  )  + scale_fill_manual("Clades",values=Clade_color) +
  xlab(substitute(paste(Delta," POC"[1]^"exp")))  +
  ylab(substitute(paste(Delta," POC"[2]^"exp")))  + scale_y_continuous(breaks = seq(-10,50,10))
# + theme(legend.position='none')
pC

jpeg(paste(path_pannel,"p4C.jpg",sep=""),width = 5200/2/resolution, height = 4000/2/resolution,res=600/2/resolution)
print(pC)
dev.off()


# Pannel D
dt_graph = data.frame(
  species = c(data1$species,data1$species),
  clade_group = c(data1$clade_group,data1$clade_group),
  category = c(rep("POC1",nrow(data1)), rep("POC2",nrow(data1))),
  value = c(data1$expressed_overused_background_POC1,data1$expressed_overused_background_POC2))

dt_graph$clade_group_facet = str_replace_all(dt_graph$clade_group," ","\n")
dt_graph$clade_group_facet = factor(dt_graph$clade_group_facet, levels = str_replace_all(names(Clade_color)," ","\n"))

pD = ggplot(dt_graph,aes(y=value,x=category,label=species,fill=clade_group))  +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = 0 ) + geom_point(size=2,pch=21,alpha=.7)+
  geom_boxplot(alpha=.1)  + facet_wrap(~clade_group_facet,scales="free_x",ncol=11)+ theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=30, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=25, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=18, family="ubuntu condensed",angle=30,vjust=0.5),
    legend.text =  element_text(color="black", size=10, family="ubuntu condensed"),
    strip.text = element_text(size=13, family="ubuntu condensed",face="bold"),
    title =  element_text(size=10, family="ubuntu condensed"),
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "grey", fill = NA)
    # strip.background = element_rect(color = "grey", size = 1)
  ) + ylab(substitute(paste(Delta," POC"^"exp"))) +
  theme(legend.position='none') + scale_fill_manual(values=Clade_color) + xlab("") 
pD

# Difference in proportion of POC between\nthe top 5% and bottom 50% expressed (%)
jpeg(paste(path_pannel,"p4D.jpg",sep=""), width = 5500/1/resolution, height = 3000/1/resolution,res=460/1/resolution)
print(pD)
dev.off()





# Figure 4

imgA = load.image(paste(path_pannel,"p4A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p4B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p4C.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p4D.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))
# Drosophila_melanogaster = readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))
# clade_png<-readPNG(paste(path_require,"clade.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure4.pdf",sep=""), width=9, height=5.6)
  
  m = matrix(rep(NA,100*100), nrow=100)
  
  for(i in 1:50){
    m[i,]=c(rep(1,48),rep(2,100-48))
  }
  
  for(i in 50:100){
    m[i,]=c(rep(3,40),rep(4,60))
  }
  layout(m)
  m
  par(mar=c(1, 1, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=0,adj=-0.2, side=2, line=1, font=2, cex=1.7,las=2)
  xhuman = 470/resolution
  yhuman = -350/resolution
  rasterImage(human,xleft=0+xhuman, ybottom=450/.9/resolution-yhuman, xright=190/.9/resolution+xhuman, ytop=0-yhuman)
  par(mar=c(0, 0, 1, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=0,adj=0.5/resolution, side=2, line=1, font=2, cex=1.7,las=2)
  xcel=500/resolution
  ycel=-350/resolution
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/1.5/resolution-ycel, xright=1000/1.5/resolution+xcel, ytop=0-ycel)
  par(mar=c(0,0, 1, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=0,adj=-2.5/resolution, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0,1, 1, 0))
  plot(imgD, axes=FALSE)
  mtext("D",at=0,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  dev.off()
}

