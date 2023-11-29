# Generate Figure 3
source("figure/figure_main_generator/library_path.R")

set_color = c( "WCp + abond" = "#33A02C" , "WCp" = "#B2DF8A" , "WBp + abond" = "#E31A1C" , "WBp" = "#FB9A99" , "not decoded" = "#e2cc1a" )

wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")


# Pannel 3 A

data4 = read.delim("data/data4.tab")
data4$Var1 = factor(data4$Var1,levels =  names(set_color))

dt_graph = data4[data4$species == "metazoa",] 

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA") 
dt_graph$codon = factor(dt_graph$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","T","A","G"),sep=""))) ) 

# dt_graph$title = factor(paste(dt_graph$codon," (",dt_graph$WB_type,")",sep=""),  
dt_graph$title = factor(paste(dt_graph$codon,sep="") , 
                        # sapply(levels(dt_graph$codon),function(x) paste(x," (",wobble_type[substr(x,3,3)],")",sep="")) )
                        sapply(levels(dt_graph$codon),function(x) paste(x,sep="")) )

dt_graph[dt_graph$amino_acid == "Ter (3)",]$Prop = NA
dt_graph[dt_graph$amino_acid == "Met (1)",]$Prop = NA
dt_graph[dt_graph$amino_acid == "Trp (1)",]$Prop = NA

pA = ggplot(dt_graph, aes(x = "" , y = Prop, fill = fct_inorder(Var1))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") + facet_wrap(~title+paste(amino_acid),nrow=4,dir="v")+
  theme_void() + theme(
    title =  element_text(size=36, family="economica"),
    legend.text =  element_text(size=40, family="economica"),
    strip.text = element_text(size=30, family="economica",face="bold"),
    legend.spacing.x = unit(1, 'cm'),
    legend.position="top",
    legend.box.spacing = unit(2, "cm"),
    plot.title = element_text(hjust = 0.5,margin = margin(0,0,20,0))
  ) +
  scale_fill_manual("",values = set_color,breaks = names(set_color))  + 
  ggtitle(paste("Metazoa"," (",dt_graph$total[1],")",sep=""))
pA

jpeg(paste(path_pannel,"p3A.jpg",sep=""), width = 8000/1,  3600/1,res=300/1)
print(pA)
dev.off()



# Pannel 3 B

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA","XX")
dt_graph = data4
dt_graph$amino_acid = str_replace_all(dt_graph$amino_acid," [(][:digit:][)]","")
dt_graph = dt_graph[dt_graph$species == "Homo_sapiens",]
dc = dt_graph[dt_graph$amino_acid == "Thr",]
dc$amino_acid = "qua"
dc$codon = c("XXA","XXT","XXC","XXG")
dt_graph = rbind(dt_graph,dc)

dt_graph$codon = factor(dt_graph$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","T","A","G"),sep=""))) )
dt_graph$title = dt_graph$codon

amino_acid_list = unique(c(dt_graph[order(dt_graph$codon)[seq(1,68,4)],]$amino_acid,dt_graph[order(dt_graph$codon)[seq(3,68,4)],]$amino_acid))

pB = ggplot(dt_graph[dt_graph$amino_acid == "Met",], aes(x = "" , y = Prop , group=amino_acid, fill = fct_inorder(Var1))) + theme_void()  
for (aa in amino_acid_list[!amino_acid_list %in% c("qua","Met","Trp","Ter")]){
  title_label = c("XXA","XXT","XXC","XXG",as.character(dt_graph[dt_graph$amino_acid == aa,]$codon))
  names(title_label) = title_label
  title_label[1:4] =  ""
  if (nrow(dt_graph[dt_graph$amino_acid == aa,]) <= 4){
    pB = pB + ggplot(dt_graph[dt_graph$amino_acid == "qua",], aes(x = "" , y = Prop , group=amino_acid, fill = fct_inorder(Var1)))
  } else {
    pB = pB + ggplot(dt_graph[dt_graph$amino_acid == aa,], aes(x = "" , y = Prop , group=amino_acid, fill = fct_inorder(Var1)))
  }
  
  pB = pB +
    ylab(aa) +
    geom_col(data=dt_graph[dt_graph$amino_acid == aa,],width = 0.1) +
    coord_polar(theta = "y") + facet_wrap(~codon,ncol=4,dir="h",labeller = labeller(codon = title_label
    ))+
    theme_void() + theme(
      title =  element_text(size=0, family="economica"),
      plot.title = element_blank(),
      legend.title =  element_text(size=20, family="economica"),
      legend.text =  element_text(size=0, family="economica"),
      strip.text = element_text(size=13, family="economica"),
      legend.position="left",
      legend.key.size = unit(0, 'cm'),
      legend.spacing.y = unit(1.5, 'cm'),
      legend.spacing.x = unit(-1, 'cm'),
      plot.margin = unit(c(0,0.3,0,0), "cm"),
      axis.title.y = element_text(size=0),
      axis.text.x =  element_text( size=0)
    ) + 
    scale_fill_manual(aa,values = set_color,breaks = names(set_color)) +
    guides(fill = guide_legend(override.aes = list(col="white",fill="white",size=0)))
}
pB

jpeg(paste(path_pannel,"p3B.jpg",sep=""), width =2000/1,  1200/1,res=270/1)
print(pB + plot_layout(ncol = 4))
dev.off()


# Pannel 3 C

dt_graph = data4
dt_graph$amino_acid = str_replace_all(dt_graph$amino_acid," [(][:digit:][)]","")
dt_graph = dt_graph[dt_graph$species == "Caenorhabditis_elegans",]
# dt_graph = dt_graph[dt_graph$species == "Drosophila_melanogaster",]
dc = dt_graph[dt_graph$amino_acid == "Thr",]
dc$amino_acid = "qua"
dc$codon = c("XXA","XXT","XXC","XXG")
dt_graph = rbind(dt_graph,dc)

dt_graph$codon = factor(dt_graph$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","T","A","G"),sep=""))) )

amino_acid_list = unique(c(dt_graph[order(dt_graph$codon)[seq(1,68,4)],]$amino_acid,dt_graph[order(dt_graph$codon)[seq(3,68,4)],]$amino_acid))


pC = ggplot(dt_graph[dt_graph$amino_acid == "Met",], aes(x = "" , y = Prop , group=amino_acid, fill = fct_inorder(Var1))) + theme_void()
for (aa in amino_acid_list[!amino_acid_list %in% c("qua","Met","Trp","Ter")]){
  title_label = c("XXA","XXT","XXC","XXG",as.character(dt_graph[dt_graph$amino_acid == aa,]$codon))
  names(title_label) = title_label
  title_label[1:4] =  ""
  if (nrow(dt_graph[dt_graph$amino_acid == aa,]) <= 4){
    pC = pC + ggplot(dt_graph[dt_graph$amino_acid == "qua",], aes(x = "" , y = Prop , group=amino_acid, fill = fct_inorder(Var1)))
  } else {
    pC = pC + ggplot(dt_graph[dt_graph$amino_acid == aa,], aes(x = "" , y = Prop , group=amino_acid, fill = fct_inorder(Var1)))
  }
  
  pC = pC +
    ylab(aa) +
    geom_col(data=dt_graph[dt_graph$amino_acid == aa,],width = 0.1) +
    coord_polar(theta = "y") + facet_wrap(~codon,ncol=4,dir="h",labeller = labeller(codon = title_label
    ))+
    theme_void() + theme(
      title =  element_text(size=0, family="economica"),
      plot.title = element_blank(),
      legend.title =  element_text(size=20, family="economica"),
      legend.text =  element_text(size=0, family="economica"),
      strip.text = element_text(size=13, family="economica"),
      legend.position="left",
      legend.key.size = unit(0, 'cm'),
      legend.spacing.y = unit(1.5, 'cm'),
      legend.spacing.x = unit(-1, 'cm'),
      plot.margin = unit(c(0,0.3,0,0), "cm"),
      axis.title.y = element_text(size=0),
      axis.text.x =  element_text( size=0)
    ) + 
    scale_fill_manual(aa,values = set_color,breaks = names(set_color)) +
    guides(fill = guide_legend(override.aes = list(col="white",fill="white",size=0)))
}
pC

jpeg(paste(path_pannel,"p3C.jpg",sep=""), width =2000/1,  1200/1,res=270/1)
print(pC + plot_layout(ncol = 4))
dev.off()


# Figure 3

imgA = load.image(paste(path_require,"wobble_pairing.png",sep="") )
imgB = load.image(paste(path_pannel,"p3A.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p3B.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p3C.jpg",sep="") )
human<-readPNG(paste(path_require,"human.png",sep=""))
Caenorhabditis_elegans = readPNG(paste(path_require,"Caenorhabditis_elegans.png",sep=""))
Drosophila_melanogaster = readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure3.pdf",sep=""), width=7, height=4.8)
  
  m=matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:5){
    m[i,]=c(rep(1,2),rep(2,8))
  }
  
  for(i in 6:10){
    m[i,]=c(rep(3,5),rep(4,5))
  }
  layout(m)
  m
  
  par(mar=c(0, 0, 1.5, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=-11,adj=-1.5, side=2, line=1, font=2, cex=1.4,las=2)
  par(mar=c(0, 0, 0.2, 0))
  plot(imgB, axes=FALSE)
  
  par(mar=c(0, 1, 1, 0))
  plot(imgC, axes=FALSE)
  mtext("B",at=-50,adj=-1, side=2, line=1, font=2, cex=1.4,las=2)
  xhuman=250
  yhuman=30
  rasterImage(human,xleft=0+xhuman, ybottom=450/2-yhuman, xright=190/2+xhuman, ytop=0-yhuman)

  plot(imgD, axes=FALSE)
  mtext("C",at=-50,adj=-1, side=2, line=1, font=2, cex=1.4,las=2)
  xcel=150
  ycel=-40
  rasterImage(Caenorhabditis_elegans,xleft=0+xcel, ybottom=350/3.5-ycel, xright=1000/3.5+xcel, ytop=0-ycel)
  dev.off()
}

