# Generate Supplementary Figure 1
source("figure/figure_supp_generator/library_path.R")
resolution = 4

# Pannel A

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
data1 = data1[ data1$nb_genes_filtered >= 5000,]

dt_graph = data1

ylabel = "var_gc3"
xlabel = "var_gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

dt_graph[,c(ylabel,xlabel)] = sqrt(dt_graph[,c(ylabel,xlabel)])

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree = arbrePhylo,display_other=F,pagels_obliged=T)

pA = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=-1, size=26,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=26, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=23, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=23, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    plot.caption = element_text(hjust = 0.370,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot",
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed",vjust = 1.5,margin = margin(l = .4,unit="cm",t=.2)),
    legend.title = element_text(color="black", size=23, family="ubuntu condensed",margin = margin(l = 0,unit="cm",t=1, b=.5)),
    legend.box.spacing =  unit(1, 'cm'),
    legend.margin =  margin(l = 0,unit="cm",t=.3)
  )  + scale_fill_manual("Clades",values=Clade_color) +
  ylab("GC3 standard deviation") +
  xlab("GCi standard deviation") +
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", p-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) +
  guides(fill = guide_legend(override.aes = list(size=5)),
  )
print(pA)

jpeg(paste(path_pannel,"p1A_supp.jpg",sep=""), width = 6000/resolution, height = 3500/resolution,res=600/resolution)
print(pA)
dev.off()


# Pannel B

ylabel = "rho_gc3_gci"
xlabel = "var_gci"

dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pB = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=-1, size=26,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=26, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=23, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=23, family="ubuntu condensed"),
    title =  element_text(color="black", size=18, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    plot.caption = element_text(hjust = .63,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot",
    legend.text =  element_text(color="black", size=20, family="ubuntu condensed",vjust = 1.5,margin = margin(l = .4,unit="cm",t=.2)),
    legend.title = element_text(color="black", size=23, family="ubuntu condensed",margin = margin(l = 0,unit="cm",t=1, b=.5)),
    legend.box.spacing =  unit(1, 'cm'),
    legend.margin =  margin(l = 0,unit="cm",t=.3)
  ) + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Spearmann rho GC3 vs GCi") +
  xlab("GCi standard deviation")  +
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", p-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) +
  guides(fill = guide_legend(override.aes = list(size=5)),
  )  + theme(legend.position='none')
print(pB)

jpeg(paste(path_pannel,"p1B_supp.jpg",sep=""), width = 4500/resolution, height = 3500/resolution,res=600/resolution)
print(pB)
dev.off()


# Supplementary Figure 1

imgA = load.image(paste(path_pannel,"p1A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p1B_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure1_supp.pdf",sep=""), width=3, height=5*0.66)
  
  m = matrix(rep(c(1,2),1*2), nrow=2)
  
  layout(m)
  m
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A", side=2,at=111/resolution, line=1, font=2, cex=1,las=2)
  par(mar=c(0, 0, 0.8, 2))
  plot(imgB, axes=FALSE)
  mtext("B", side=2,adj=-1.6,at=0, line=1, font=2, cex=1,las=2)
  dev.off()
}
