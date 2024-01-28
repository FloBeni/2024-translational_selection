# Generate Figure 8
source("figure/figure_main_generator/library_path.R")


# Pannel 8 A
data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1 = data1[data1$clade_group %in% c("Diptera","Lepidoptera"),]

dt_graph = data1

ylabel = "g_abond_ag"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
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
    caption = substitute(paste(model," :",aic," R"^2,"= ",r2,", p-value = ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  )  + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Average GC of most abundant tRNAs set\n(Unn/Cnn)") + 
  xlab("Average per gene GCi") 

pA

jpeg(paste(path_pannel,"p8A.jpg",sep=""),width = 5200/2, height = 4000/2,res=700/2)
print(pA)
dev.off()


# Pannel 7 B
dt_graph = data1

ylabel = "c_duc_ic"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

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
  # geom_boxplot(data=dt_graph , aes_string(x="group", y= ylabel,group="group"),outlier.shape = NA ,pch=21,fill=set_color[5],size=1)+
  # geom_point(data=average_data , aes(x=x, y= y, group="moyenne"  ),pch=21,fill=set_color[7],size=8)+
  # geom_errorbarh(data=average_data ,aes( xmin=x-sdx,xmax=x+sdx,x=x,y=y,group="moyenne"),height=0.02, size=1,col=set_color[7])+
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1,margin = margin(t = 5)),
    legend.title = element_text(color="black", size=20, family="economica"),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + 
  labs(
    caption = substitute(paste(model," :",aic," R"^2,"= ",r2,", p-value = ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) + theme(legend.position='none') + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Average GC of most selected codons set\n(nnC/nnU)") + 
  xlab("Average per gene GCi") +   theme(legend.position = c(0.85, 0.2),
                                         legend.background = element_rect(fill="NA"),
                                         legend.spacing.x = unit(0.1, 'cm'),
                                         legend.spacing.y = unit(0.1, 'cm'),
                                         legend.box.background = element_rect(colour = "black")
  )

pB

jpeg(paste(path_pannel,"p8B.jpg",sep=""),width = 5200/2, height = 4000/2,res=700/2)
print(pB)
dev.off()


# Figure 7

imgA = load.image(paste(path_pannel,"p8A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p8B.jpg",sep="") )
imgC = load.image(paste(path_require,"hypothesis.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure8.pdf",sep=""), width=10, height=8*2/2)
  
  m = matrix(rep(NA,2*2), nrow=2)
  
  m[1,]=c(1,2)
  m[2,]=c(3,3)
  layout(m)
  
  par(mar=c(0,2, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=-100,adj=0, side=2, line=1, font=2, cex=2,las=2)
  xaxis=2230
  yaxis=150
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/5+xaxis, ytop=-900/5+yaxis)
  
  plot(imgB, axes=FALSE)
  mtext("B",at=-100,adj=0, side=2, line=1, font=2, cex=2,las=2)
  xaxis=2230
  yaxis=150
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/5+xaxis, ytop=-900/5+yaxis)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgC, axes=FALSE)
  mtext("C",at=100,adj=-1, side=2, line=1, font=2, cex=2,las=2)
  
  dev.off()
}
