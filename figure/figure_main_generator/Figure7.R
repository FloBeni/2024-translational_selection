# Generate Figure 7
source("figure/figure_main_generator/library_path.R")


# Pannel 7 A
data1 = read.delim("data/data1.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1 = data1[data1$clade_group == "Mecopterida",]

dt_graph = data1

ylabel = "gc_abond_ag"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species = dt_graph$species,
                                         pgls_x = lm_x,
                                         pgls_y = lm_y), species, vcv=TRUE)


pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn), list(nbspecies=nrow(dt_graph),
                                                           lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                           pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  )  + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Average GC of most abundant tRNAs A/G") + 
  xlab("Average per gene GCi")

pA

jpeg(paste(path_pannel,"F7pA.jpg",sep=""),width = 5200/2, height = 4000/2,res=700/2)
print(pA)
dev.off()


# Pannel 7 B
dt_graph = data1

ylabel = "gc_duc"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)

average_data = data.frame(x=tapply(dt_graph[,xlabel],dt_graph[,ylabel],mean),
                          y=tapply(dt_graph[,ylabel],dt_graph[,ylabel],mean),
                          sdx=tapply(dt_graph[,xlabel],dt_graph[,ylabel],sd),
                          sdy=tapply(dt_graph[,ylabel],dt_graph[,ylabel],sd)
)

dt_graph[,"group"] = tapply(dt_graph[,xlabel],dt_graph[,xlabel]<=median(dt_graph[,xlabel]) , mean)[as.character(dt_graph[,xlabel]<=median(dt_graph[,xlabel]))]

x=dt_graph[dt_graph$group == unique(dt_graph$group)[1],][,ylabel]
y=dt_graph[dt_graph$group == unique(dt_graph$group)[2],][,ylabel]
wilcox.test(x,y)
t.test(x,y)
t.test(x, y, var.equal = FALSE)


pB =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_boxplot(data=dt_graph , aes_string(x="group", y= ylabel,group="group"),outlier.shape = NA ,pch=21,fill=set_color[5],size=1)+
  # geom_point(data=average_data , aes(x=x, y= y, group="moyenne"  ),pch=21,fill=set_color[7],size=8)+
  # geom_errorbarh(data=average_data ,aes( xmin=x-sdx,xmax=x+sdx,x=x,y=y,group="moyenne"),height=0.02, size=1,col=set_color[7])+
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn), list(nbspecies=nrow(dt_graph),
                                                           lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                           pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Average GC of most selected codon set") + 
  xlab("Average per gene GCi")

pB

jpeg(paste(path_pannel,"F7pB.jpg",sep=""),width = 5200/2, height = 4000/2,res=700/2)
print(pB)
dev.off()


# Figure 7

imgA = load.image(paste(path_pannel,"F7pA.jpg",sep="") )
imgB = load.image(paste(path_pannel,"F7pB.jpg",sep="") )
imgC = load.image(paste(path_require,"hypothesis.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure7.pdf",sep=""), width=10, height=8)
  
  m = matrix(rep(NA,2*2), nrow=2)
  
  m[1,]=c(1,2)
  m[2,]=c(3,3)
  layout(m)
  
  par(mar=c(0,2, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=-100,adj=0, side=2, line=1, font=2, cex=2,las=2)
  xaxis=2230
  yaxis=1430
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/4.5+xaxis, ytop=-900/4.5+yaxis)
  
  plot(imgB, axes=FALSE)
  mtext("B",at=-100,adj=0, side=2, line=1, font=2, cex=2,las=2)
  xaxis=2230
  yaxis=1430
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/4.5+xaxis, ytop=-900/4.5+yaxis)
  
  plot(imgC, axes=FALSE)
  mtext("C",at=10,adj=-5, side=2, line=1, font=2, cex=2,las=2)
  dev.off()
}
