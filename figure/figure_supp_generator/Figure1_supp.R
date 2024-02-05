# Generate Supplementary Figure 1
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 1 A

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
data1 = data1[ data1$nb_genes_filtered >= 5000,]

dt_graph = data1

ylabel = "var_gc3"
xlabel = "var_gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

dt_graph[,c(ylabel,xlabel)] = sqrt(dt_graph[,c(ylabel,xlabel)])

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pA = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 5)),
    legend.title =  element_text(color="black", size=25, family="economica"),
    plot.caption = element_text(hjust = 0.370, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )  + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Standard deviation per gene GC3") +
  xlab("Standard deviation per gene GCi") +
  labs(
    caption = substitute(paste(model," :",aic," R"^2,"= ",r2,", p-value = ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) +
  guides(fill = guide_legend(override.aes = list(size=5)),
  )
print(pA)

jpeg(paste(path_pannel,"p1A_supp.jpg",sep=""), width = 6000/1, height = 3500/1,res=600/1)
print(pA)
dev.off()


# Supplementary Pannel 1 B

ylabel = "rho_gc3_gci"
xlabel = "var_gci"

dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pB = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  # geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 7)),
    plot.caption = element_text(hjust = 0.6, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  ) + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Spearmann Rho GC3 GCi per species") +
  xlab("Standard deviation per gene GCi")  +
  labs(
    caption = substitute(paste(model," :",aic," R"^2,"= ",r2,", p-value = ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) +
  guides(fill = guide_legend(override.aes = list(size=5)),
  )  + theme(legend.position='none')
print(pB)

jpeg(paste(path_pannel,"p1B_supp.jpg",sep=""), width = 4500/1, height = 3500/1,res=600/1)
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
  mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  par(mar=c(0, 0, 0.8, 2))
  plot(imgB, axes=FALSE)
  mtext("B", side=2,adj=-1.6,at=0, line=1, font=2, cex=1,las=2)
  dev.off()
}
