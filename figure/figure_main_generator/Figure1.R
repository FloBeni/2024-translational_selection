source("figure/figure_main generator/library_path.R")

############## Pannel 1 A
arbrePhylotips = arbrePhylo
arbrePhylotips$tip.label <- str_replace_all(arbrePhylotips$tip.label,"_"," ")
edge_group <- str_replace_all(arbrePhylotips$tip.label,"_"," ")
edge_clade <- rep("branch",length(arbrePhylotips$edge[,2]))
for (group in unique(edge_group)){
  if (group %in% unlist(listNomSpecies)){
    edge_clade[arbrePhylotips$edge[,2] %in% grep(group,edge_group)] =
      names(listNomSpecies[unlist(lapply(listNomSpecies,function(x) group %in% x))])
  }
}

for (clade in  names(listNomSpecies)){print(clade)
  edge_clade[ which.edge(arbrePhylotips, arbrePhylotips$edge[,2][edge_clade == clade] ) ] = clade
}
node_metadata = data.frame(node=arbrePhylotips$edge[,2],color=edge_clade)
node_metadata$color = factor(node_metadata$color, levels = c("Mecopterida","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Teleostei","Mammalia","Aves","Other Tetrapods"))
p1A = ggtree(arbrePhylotips,layout="roundrect",size=1)
# %>% flip(264, 375)
p1A <- p1A %<+% node_metadata  + aes(color=color) + 
  scale_color_manual("Clade",values=Clade_color[unique(edge_clade)]) + theme(
    panel.background = element_rect(fill = "white", linetype = "dashed")
  ) + theme(legend.position = "none")
p1A

# jpeg(paste(path_pannel,"p1A.jpg",sep=""), width = 1500/1, height = 3000/1,res=250/1)
# print(p1A)
# dev.off()


############## Pannel 1 B
data1 = read.delim("data/data1_bis.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

dt_graph = data1

ylabel = "gc3"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)

p1B = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  # geom_errorbar(aes(ymin = gc3-sqrt(var_gc3), ymax = gc3+sqrt(var_gc3)),alpha=0.2) +
  # geom_errorbarh(aes(xmin = gci-sqrt(var_gci), xmax =gci+sqrt(var_gci)),alpha=0.2) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.65, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Average per gene GC3") +
  xlab("Average per gene GCi") +
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn," / PGLS:"," R"^2,pgls_eq), list(nbspecies=nrow(dt_graph),
                                                                                     lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                                                     pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird)))),
    title = substitute(paste("N = ",nbspecies," species"), list(nbspecies=nrow(dt_graph),
                                                            lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                            pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  ) 
print(p1B)



jpeg(paste(path_pannel,"p1B.jpg",sep=""), width = 4000/1, height = 3000/1,res=500/1)
print(p1B)
dev.off()

############## Pannel 1 C
data2 = read.delim("data/data2_bis.tab")
dt_graph = data2[data2$species == "Homo_sapiens" ,]
spearman_method_aa = cor.test( dt_graph$GCi, dt_graph$GC3,method="spearman",exact=F)


p1C = ggplot(dt_graph ,
             aes(x=GCi ,y=GC3))  +
  geom_point(col=set_color[8],alpha=0.4,size=.51)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(24,22,21,23,25,20))+ xlab("log10( FPKM+1 )") + ylab("Frequency optimal codons") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=30,family="economica"),
    axis.title.y = element_text(color="black", size=30, family="economica"),
    axis.text.y =  element_text(color="black", size=25, family="economica"),
    axis.text.x =  element_text(color="black", size=25, family="economica"),
    title =  element_text(color="black", size=18, family="economica"),
    legend.text =  element_text(color="black", size=16, family="economica"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.65, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )  +labs(fill='Categories',color='Categories',shape='',linetype='')+ 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = spearman_method_aa$p.value))
  ) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2,size=NA),
  )  +  ylab("GC3 rate per gene") + xlab("GCi rate per gene") +xlim(0.1,.8) +ylim(0.15,1)
p1C

jpeg(paste(path_pannel,"p1C.jpg",sep=""), 
     width = 6000/4.5,  5000/4,res=1000/4)
print(p1C)
dev.off()


############## Pannel 3 D

dt_graph = data2[data2$species == "Caenorhabditis_elegans" ,]
spearman_method_aa = cor.test( dt_graph$GCi, dt_graph$GC3,method="spearman",exact=F)


p1D = ggplot(dt_graph ,
             aes(x=GCi ,y=GC3))  +
  geom_point(col=set_color[8],alpha=0.4,size=.51)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(24,22,21,23,25,20)) + xlab("log10( FPKM+1 )") + ylab("Frequency optimal codons") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=30,family="economica"),
    axis.title.y = element_text(color="black", size=0, family="economica"),
    axis.text.y =  element_text(color="black", size=0, family="economica"),
    axis.text.x =  element_text(color="black", size=25, family="economica"),
    title =  element_text(color="black", size=18, family="economica"),
    legend.text =  element_text(color="black", size=16, family="economica"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.55, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+labs(fill='Categories',color='Categories',shape='',linetype='') + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value = ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  ) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2,size=NA),
  )  +  ylab("") + xlab("GCi rate per gene") +xlim(0.1,.8) +ylim(0.15,1)
p1D

jpeg(paste(path_pannel,"p1D.jpg",sep=""), 
     width = 5000/4.5,  5000/4,res=1000/4)
print(p1D)
dev.off()



############## Figure 1 

imgA = load.image(paste(path_pannel,"p1A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p1B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p1C.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p1D.jpg",sep="") )
clade_png<-readPNG(paste(path_require,"clade.png",sep=""))
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))


aves<-readPNG(paste(path_require,"aves.png",sep=""))
teleostei<-readPNG(paste(path_require,"teleostei.png",sep=""))
monkey<-readPNG(paste(path_require,"monkey.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))
bee<-readPNG(paste(path_require,"bee.png",sep=""))


{
  pdf(file= paste(path_figure,"Figure1.pdf",sep=""), width=6.5, height=5)
  
  m=matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:10){
    m[,i]=rep(1)
  }
  
  for(i in 5:10){
    m[,i]=c(rep(2,5),rep(3,5))
  }
  
  for(i in 6:10){
    m[i,]=c(rep(1,4),rep(3,3),rep(4,3))
  }
  layout(m)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgA, axes=FALSE)
  mtext("A",at=49.4,adj=-1, side=2, line=1, font=2, cex=1.4,las=2)
  
  xmonkey=920
  ymonkey=2000
  rasterImage(clade_png,xleft=0+xmonkey, ybottom=800/0.85+ymonkey, xright=500/.85+xmonkey, ytop=ymonkey)
  
  xaxis=700
  yaxis=2800
  rasterImage(teleostei,xleft=0+xaxis, ybottom=0+yaxis, xright=900/6+xaxis, ytop=-500/6+yaxis)

  xaxis=680
  yaxis=2350
  rasterImage(aves,xleft=0+xaxis, ybottom=0+yaxis, xright=600/3.5+xaxis, ytop=-750/3.5+yaxis)
  
  xaxis=750
  yaxis=1750
  rasterImage(monkey,xleft=0+xaxis, ybottom=0+yaxis, xright=900/5+xaxis, ytop=-900/5+yaxis)
  
  xcel=1300
  ycel=1050
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=0+ycel, xright=1000/5+xcel, ytop=-350/5+ycel)
  
  xaxis=1100
  yaxis=750
  rasterImage(bee,xleft=0+xaxis, ybottom=0+yaxis, xright=900/6+xaxis, ytop=-700/6+yaxis)
  
  xaxis=1200
  yaxis=350
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=900/8+xaxis, ytop=-900/8+yaxis)
  
  par(mar=c(0, 1, 2, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=49.4,adj=-1.5, side=2, line=1, font=2, cex=1.4,las=2)
  
  par(mar=c(0, 0, 0, 0))
  plot(imgC, axes=FALSE)
  xhuman=300
  yhuman=-60
  rasterImage(human,xleft=0+xhuman, ybottom=350/1.7-yhuman, xright=190/1.7+xhuman, ytop=0-yhuman)
  mtext("C", adj=0.1,side=2,at=0, line=1, font=2, cex=1.4,las=2)
  par(mar=c(0, 0, 0, 2.5))
  plot(imgD, axes=FALSE)
  xcel=125
  ycel=-50
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/3-ycel, xright=1000/3+xcel, ytop=0-ycel)
  dev.off()
  
}

