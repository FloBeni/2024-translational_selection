# Generate Figure 6
source("figure/figure_main_generator/library_path.R")
resolution = 1

# Pannel A
data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1 = data1[data1$clade_group %in% c("Diptera","Lepidoptera"),]
data1 = data1[data1$species != "Eumeta_japonica",]

dt_graph = data1

ylabel = "g_abond_ag"
xlabel = "gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() +  theme(
    axis.title.x = element_text(color="black", size=25,vjust=0,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=25,vjust=1.5, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=23, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=23, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed",vjust = 1,margin = margin(t = 5)),
    legend.title = element_text(color="black", size=20, family="ubuntu condensed"),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", p-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  )  + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Average CNN of preferred\nisodecoder tRNAs (UNN/CNN)") + 
  xlab("GCi") 

pA

jpeg(paste(path_pannel,"p6A.jpg",sep=""),width = 5200/2/resolution, height = 4000/2/resolution,res=700/2/resolution)
print(pA)
dev.off()


# Pannel B
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
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,vjust=0,family="ubuntu condensed"),
    axis.title.y = element_text(color="black", size=25,vjust=1.5, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=23, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=23, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed",vjust = 1,margin = margin(t = 5)),
    legend.title = element_text(color="black", size=20, family="ubuntu condensed"),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + 
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", p-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) + theme(legend.position='none') + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Average NNC of\npreferred codons (NNT/NNC)") + 
  xlab("GCi") +   theme(legend.position = c(0.84, 0.15),
                                         legend.background = element_rect(fill="NA"),
                                         legend.spacing.x = unit(0.1, 'cm'),
                                         legend.spacing.y = unit(0.1, 'cm'),
                                         legend.box.background = element_rect(colour = "black")
  )

pB

jpeg(paste(path_pannel,"p6B.jpg",sep=""),width = 5200/2/resolution, height = 4000/2/resolution,res=700/2/resolution)
print(pB)
dev.off()


# Figure 6

imgA = load.image(paste(path_pannel,"p6A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p6B.jpg",sep="") )
imgC = load.image(paste(path_require,"hypothesis.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure6.pdf",sep=""), width=10, height=8*2/2)
  
  m = matrix(rep(NA,2*2), nrow=2)
  
  m[1,]=c(1,2)
  m[2,]=c(3,3)
  layout(m)
  
  par(mar=c(0,2, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=-100,adj=0, side=2, line=1, font=2, cex=2.3,las=2)
  xaxis=2230
  yaxis=150
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/5+xaxis, ytop=-900/5+yaxis)
  
  plot(imgB, axes=FALSE)
  mtext("B",at=-100,adj=0, side=2, line=1, font=2, cex=2.3,las=2)
  xaxis=2230
  yaxis=150
  rasterImage(fly,xleft=0+xaxis, ybottom=0+yaxis, xright=1200/5+xaxis, ytop=-900/5+yaxis)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgC, axes=FALSE)
  mtext("C",at=100,adj=-1, side=2, line=1, font=2, cex=2.3,las=2)
  
  dev.off()
}
