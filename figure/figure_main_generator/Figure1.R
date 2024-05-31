# Generate Figure 1
source("figure/figure_main_generator/library_path.R")

# Pannel A

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

edge_clade_prev = edge_clade
list_inclusion =  list("Other Metazoans"=c("Diptera","Lepidoptera","Coleoptera","Other Insects","Other Tetrapods","Nematoda","Hymenoptera","Teleostei","Other Metazoans"),
                       "Other Tetrapods"=c("Aves","Mammalia","Other Tetrapods"),"Other Insects"=c("Diptera","Lepidoptera","Coleoptera","Other Insects","Hymenoptera"),
                       Nematoda="Nematoda",Teleostei="Teleostei",Hymenoptera="Hymenoptera",Aves="Aves",Mammalia="Mammalia","Diptera"="Diptera","Lepidoptera"="Lepidoptera","Coleoptera"="Coleoptera"
)


clade="Other Metazoans"

for (clade in names(list_inclusion)){
  edge_clade[ which.edge(arbrePhylotips,  arbrePhylotips$edge[,2][edge_clade_prev %in% unlist(list_inclusion[clade])] ) ] = clade
}
node_metadata = data.frame(node=arbrePhylotips$edge[,2],color=edge_clade)

node_metadata$color = factor(node_metadata$color, levels = names(Clade_color))

offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree") # to use flip deprecated

pA = ggtree(arbrePhylotips,layout="roundrect",size=1) %>% flip(264, 375)
pA <- pA %<+% node_metadata  + aes(color=color) + 
  scale_color_manual("Clade",values=Clade_color[unique(edge_clade)]) + theme(
    panel.background = element_rect(fill = "white", linetype = "dashed")
  ) + theme(legend.position = "none")

pA

jpeg(paste(path_pannel,"p1A.jpg",sep=""), width = 1500/1, height = 3000/1,res=250/1)
print(pA)
dev.off()


# Pannel B

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

labels = paste(names(tapply(data1$species,data1$clade_group,length))," N=",tapply(data1$species,data1$clade_group,length),sep="")
names(labels) = names(tapply(data1$species,data1$clade_group,length))

data1 = data1[ data1$nb_genes_filtered >= 5000,]
dt_graph = data1

ylabel = "gc3"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pB = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  geom_abline(linetype="dashed") +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=28,vjust=0.5,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=28,vjust=1.5, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=24, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=24, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=24, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )  + scale_fill_manual("Clades",values=Clade_color,labels=labels) +
  ylab("GC3") +
  xlab("GCi") +
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", p-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  )+ theme(legend.position='none')
pB



jpeg(paste(path_pannel,"p1B.jpg",sep=""), width = 4000/3, height = 3500/3,res=600/3)
print(pB)
dev.off()


# Pannel C

data2 = read.delim("data/data2_supp.tab")
dt_graph = data2[data2$species == "Homo_sapiens" ,]
spearman_method_aa = cor.test( dt_graph$GCi, dt_graph$GC3,method="spearman",exact=F)
if(spearman_method_aa$p.value < 1e-16){
  spearman_method_aa$p.value = "< 1e-16"
} else {
  spearman_method_aa$p.value = paste("= ",spearman_method_aa$p.value,sep="")
}

pC = ggplot(dt_graph , aes(x=GCi ,y=GC3))  +
  geom_abline(linetype="dashed") +
  geom_point(col=set_color[8],alpha=0.4,size=.51)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(24,22,21,23,25,20))+ xlab("log10( FPKM+1 )") + ylab("Frequency optimal codons") + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=0.5, size=30,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=30,vjust=1.5, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=25, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=25, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=16, family="ubuntu condensed"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.7, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )  +labs(fill='Categories',color='Categories',shape='',linetype='')+ 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = spearman_method_aa$p.value))
  ) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2,size=NA),
  )  +  ylab("GC3") + xlab("GCi") + xlim(0.1,.8) +ylim(0.15,1)
pC

jpeg(paste(path_pannel,"p1C.jpg",sep=""), 
     width = 6000/4.5,  5000/4,res=1000/4)
print(pC)
dev.off()


# Pannel D

dt_graph = data2[data2$species == "Caenorhabditis_elegans" ,]
spearman_method_aa = cor.test( dt_graph$GCi, dt_graph$GC3,method="spearman",exact=F)
if(spearman_method_aa$p.value < 1e-16){
  spearman_method_aa$p.value = "< 1e-16"
} else {
  spearman_method_aa$p.value = paste("= ",spearman_method_aa$p.value,sep="")
}



pD = ggplot(dt_graph , aes(x=GCi ,y=GC3))  +
  geom_abline(linetype="dashed") +
  geom_point(col=set_color[8],alpha=0.4,size=.51)+
  scale_fill_manual(values=set_color) +
  scale_color_manual(values=set_color) +
  scale_shape_manual(values=c(24,22,21,23,25,20)) + xlab("log10( FPKM+1 )") + ylab("Frequency optimal codons") + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=0.5, size=30,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=0, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=0, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=25, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=16, family="ubuntu condensed"),
    strip.text = element_text(size=15),
    plot.caption = element_text(hjust = 0.55, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ labs(fill='Categories',color='Categories',shape='',linetype='') + 
  labs(
    caption = substitute(paste("rho = ",rho_aa_fpkm,", p-value ",pval_aa_fpkm), list(
      rho_aa_fpkm = round(spearman_method_aa$estimate, 2),
      pval_aa_fpkm = formatC(spearman_method_aa$p.value, format = "e", digits = 0)))
  ) +
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2,size=NA),
  )  +  ylab("") + xlab("GCi") + xlim(0.1,.8) +ylim(0.15,1)
pD

jpeg(paste(path_pannel,"p1D.jpg",sep=""), 
     width = 5000/4.5,  5000/4,res=1000/4)
print(pD)
dev.off()



# Figure 1 

imgA = load.image(paste(path_pannel,"p1A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p1B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p1C.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p1D.jpg",sep="") )
clade_png<-readPNG(paste(path_require,"clade.png",sep=""))
human<-readPNG(paste(path_require,"human_f.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))


aves<-readPNG(paste(path_require,"aves.png",sep=""))
teleostei<-readPNG(paste(path_require,"teleostei.png",sep=""))
monkey<-readPNG(paste(path_require,"monkey.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))
bee<-readPNG(paste(path_require,"bee.png",sep=""))
coleoptera<-readPNG(paste(path_require,"coleoptera.png",sep=""))
lepidoptera<-readPNG(paste(path_require,"lepidoptera.png",sep=""))



{
  pdf(file= paste(path_figure,"Figure1.pdf",sep=""), width=7, height=5)
  
  m=matrix(rep(NA,10*100), nrow=10)
  
  for(i in 1:100){
    m[,i]=rep(1)
  }
  
  for(i in 50:100){
    m[,i]=c(rep(2,5),rep(3,5))
  }
  
  for(i in 6:10){
    m[i,]=c(rep(1,42),rep(3,29),rep(4,29))
  }
  
  
  layout(m)
  m
  par(mar=c(0, 1, 0, 4))
  plot(imgA, axes=FALSE)
  mtext("A",at=100,adj=-1, side=2, line=1, font=2, cex=1.4,las=2)
  
  xclade=900
  yclade=1700
  rasterImage(clade_png,xleft=0+xclade, ybottom=850/.75+yclade, xright=520/.75+xclade, ytop=yclade)

  xaxis=700
  yaxis=2800
  rasterImage(teleostei,xleft=0+xaxis, ybottom=0+yaxis, xright=1100/6+xaxis, ytop=-500/6+yaxis)

  xaxis=680
  yaxis=2350
  rasterImage(aves,xleft=0+xaxis, ybottom=0+yaxis, xright=600/3.5+xaxis, ytop=-750/3.5+yaxis)

  xaxis=750
  yaxis=1750
  rasterImage(monkey,xleft=0+xaxis, ybottom=0+yaxis, xright=900/5+xaxis, ytop=-900/5+yaxis)

  xcel=1310
  ycel=1120
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=0+ycel, xright=1000/5+xcel, ytop=-350/5+ycel)

  xaxis=1100
  yaxis=800
  rasterImage(bee,xleft=0+xaxis, ybottom=0+yaxis, xright=900/5+xaxis, ytop=-700/5+yaxis)

  xaxis=1150
  yaxis=440
  rasterImage(coleoptera,xleft=0+xaxis, ybottom=0+yaxis, xright=1500/10+xaxis, ytop=-900/10+yaxis)

  xaxis=1250
  yaxis=350
  rasterImage(lepidoptera,xleft=0+xaxis, ybottom=0+yaxis, xright=1500/7+xaxis, ytop=-900/7+yaxis)

  xaxis=1300
  yaxis=150
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/10+xaxis, ytop=-900/10+yaxis)

  par(mar=c(0, 0, 1, 2.5))
  plot(imgB, axes=FALSE)
  mtext("B",at=49.4,adj=-1.5, side=2, line=1, font=2, cex=1.4,las=2)

  par(mar=c(0, 0, 0, 0))
  plot(imgC, axes=FALSE)
  xhuman=300
  yhuman=-60
  rasterImage(human,xleft=0+xhuman, ybottom=500/1.9-yhuman, xright=190/1.9+xhuman, ytop=0-yhuman)
  mtext("C", adj=0.1,side=2,at=0, line=1, font=2, cex=1.4,las=2)
  par(mar=c(0, 0, 0, 2.5))
  plot(imgD, axes=FALSE)
  xcel=125
  ycel=-50
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/3-ycel, xright=1000/3+xcel, ytop=0-ycel)
  dev.off()
}

