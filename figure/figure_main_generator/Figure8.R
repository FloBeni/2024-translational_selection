# Generate Figure 8
source("figure/figure_main_generator/library_path.R")


# Pannel 8 A

data1 = read.delim("data/data1.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1[,c("lifespan_days","length_cm","weight_kg")] = GTDrift_list_species[data1$species,c("lifespan_days","length_cm","weight_kg")]

dt_graph = data1
ylabel = "expressed_overused_background_global"
xlabel = "lifespan_days"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = log10(dt_graph[,ylabel])
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)


pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn," / PGLS:"," R"^2,pgls_eq), list(nbspecies=nrow(dt_graph),
                                                                                     lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                                                     pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird)))),
    title = substitute(paste("N = ",nbspecies," species"), list(nbspecies=nrow(dt_graph),
                                                            lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                            pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  ) + scale_fill_manual(values=Clade_color) +
  ylab("Translational selection intensity") + 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000),labels=c(0.05,0.1,0.5,1,5,10,100,1000,10000)) + 
  xlab("Longevity (days log scale)")+ annotation_logticks(sides="b") 

pA

jpeg(paste(path_pannel,"F8pA.jpg",sep=""),width = 4200/2, height = 4000/2,res=700/2)
print(pA)
dev.off()


# Pannel 8 B

dt_graph = data1

ylabel = "expressed_overused_background_global"
xlabel = "var_gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)


pB =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.59, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + 
  labs(
    title = substitute(paste("N = ",nbspecies," species"), list(nbspecies=nrow(dt_graph),
                                                            lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                            pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  )+ scale_fill_manual("Clades",values=Clade_color) +
  ylab("Translational selection intensity") + 
  xlab("Variance per gene GCi")

pB

jpeg(paste(path_pannel,"F8pB.jpg",sep=""),width = 6200/2, height = 3550/2,res=700/2)
print(pB)
dev.off()


# Figure 8

imgA = load.image(paste(path_pannel,"F8pA.jpg",sep="") )
imgB = load.image(paste(path_pannel,"F8pB.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure8.pdf",sep=""), width=10, height=4)
  
  m = matrix(rep(NA,1*100), nrow=1)
  
  m[1,]=c(rep(1,40),rep(2,60))
  layout(m)
  
  par(mar=c(0,2, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=0,adj=0, side=2, line=1, font=2, cex=2,las=2)
  
  par(mar=c(2,2, 2, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=0,adj=0, side=2, line=1, font=2, cex=2,las=2)
  dev.off()
}
