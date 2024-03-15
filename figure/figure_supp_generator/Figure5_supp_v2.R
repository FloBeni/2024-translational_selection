# Generate Supplementary Figure 3
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 3 A

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0 & data1$pval_aa_fpkm < 0.05  & data1$nb_genes_filtered > 5000,]

dt_graph = data1
ylabel = "expressed_overused_background_POCs"
xlabel = "expressed_overused_POCs"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 


model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_abline(linetype="dashed") +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = 0 ) +
  geom_vline(size=1,linetype="dashed",col="red", xintercept = 0 ) +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica",margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    legend.title =  element_text(color="black", size=25, family="economica"),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5)))+
  labs(
    caption = substitute(paste(model," :",aic," R"^2,"= ",r2,", p-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) + scale_fill_manual("Clades",values=Clade_color) + xlab(substitute(paste(Delta," POC"^"exp",", not accounting for intron control")))  + ylab(substitute(paste(Delta," POC"^"exp")))


pA


jpeg(paste(path_pannel,"p3A_supp.jpg",sep=""), width = 4500/1, height = 3000/1,res=400/1)
print(pA)
dev.off()


# Supplementary Figure 3

imgA = load.image(paste(path_pannel,"p3A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure5_supp.pdf",sep=""), width=6, height=4)
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}
