# Generate Supplementary Figure 3
source("figure/figure_supp_generator/library_path.R")


# Supplementary Pannel 3 A

data1 = read.delim("data/data1_supp.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0 & data1$pval_aa_fpkm < 0.05  & data1$nb_genes_filtered > 5000,]

dt_graph = data1
ylabel = "expressed_overused_WB_WC_notambiguous"
xlabel = "expressed_overused_background_WB_WC_notambiguous"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)


pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel))  +
  geom_point(aes(fill=clade_group),size=4,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26, family="economica",margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.45, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) +
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn," / PGLS:"," R"^2,pgls_eq), list(nbspecies=nrow(dt_graph),
                                                                                     lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                                                     pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird)))),
    title = substitute(paste("N = ",nbspecies," species"), list(nbspecies=nrow(dt_graph),
                                                            lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                            pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  ) + scale_fill_manual("Clades",values=Clade_color) + xlab("Difference in POC proportion between the top 5% and bottom 50% expressed (%)") +
  ylab("Difference in POC proportion between the top 5% and\nbottom 50% expressed (%, accounting for POCMT variations)")

pA


jpeg(paste(path_pannel,"p3A_supp.jpg",sep=""), width = 5500/1, height = 3000/1,res=400/1)
print(pA)
dev.off()


# Supplementary Figure 3

imgA = load.image(paste(path_pannel,"p3A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure3_supp.pdf",sep=""), width=6.75, height=4)
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}
