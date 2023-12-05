# Generate Supplementary Figure 10
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 10 A

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
dnds = read.delim("data/GTDrift_Metazoa_dNdS.tab")
rownames(dnds) = dnds$species

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1[,c("lifespan_days","length_cm","weight_kg")] = GTDrift_list_species[data1$species,c("lifespan_days","length_cm","weight_kg")]
data1[,c("dNdS")] = dnds[data1$species,c("dNdS")]

dt_graph = data1
ylabel = "expressed_overused_background_global"
xlabel = "length_cm"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = log10(dt_graph[,xlabel])
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)


pA = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.99, face= "italic", size=20, family="economica"),
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
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000,5000),labels=c(0.01,0.1,1,10,100,1000,5000)) + xlab("Body length (cm, log scale)")+ annotation_logticks(sides="b") 
pA

jpeg(paste(path_pannel,"F10pA_supp.jpg",sep=""),width = 4200/2, height = 4000/2,res=700/2)
print(pA)
dev.off()


# Supplementary Pannel 10 B

dt_graph = data1

ylabel = "expressed_overused_background_global"
xlabel = "weight_kg"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = log10(dt_graph[,xlabel])
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
    plot.caption = element_text(hjust = 0.15, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + 
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn," / PGLS:"," R"^2,pgls_eq), list(nbspecies=nrow(dt_graph),
                                                                                     lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                                                     pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird)))),
    title = substitute(paste("N = ",nbspecies," species"), list(nbspecies=nrow(dt_graph),
                                                                lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                                pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  ) + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Translational selection intensity") +  
  scale_x_log10(breaks=c(10^-6,10^-4,10^-2,10^0,10^2,10^4,10^6),labels=label_log(digits = 2),limits = c(0.000001,1000000)) + xlab("Body Weight (kg, log scale)")
pB

jpeg(paste(path_pannel,"F10pB_supp.jpg",sep=""),width = 6200/2, height = 3950/2,res=700/2)
print(pB)
dev.off()


# Supplementary Pannel 10 C

dt_graph = data1

ylabel = "expressed_overused_background_global"
xlabel = "dNdS"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)


pC = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica"),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.99, face= "italic", size=20, family="economica"),
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
  ylab("Translational selection intensity") + xlab("dN/dS")
pC

jpeg(paste(path_pannel,"F10pC_supp.jpg",sep=""),width = 4200/2, height = 4000/2,res=700/2)
print(pC)
dev.off()


# Supplementary Figure 10

imgA = load.image(paste(path_pannel,"F10pA_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"F10pB_supp.jpg",sep="") )
imgC = load.image(paste(path_pannel,"F10pC_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure10_supp.pdf",sep=""), width=10, height=8)
  
  m = matrix(rep(NA,2*100), nrow=2)
  
  m[1,]=c(rep(1,40),rep(2,60))
  m[2,]=c(rep(3,40),rep(4,60))
  layout(m)
  
  par(mar=c(0,2, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=0,adj=0, side=2, line=1, font=2, cex=2,las=2)
  
  plot(imgB, axes=FALSE)
  mtext("B",at=0,adj=0, side=2, line=1, font=2, cex=2,las=2)
  
  plot(imgC, axes=FALSE)
  mtext("C",at=0,adj=0, side=2, line=1, font=2, cex=2,las=2)
  dev.off()
}
