# Generate Figure 6
source("figure/figure_main_generator/library_path.R")


# Pannel 6 A

data1 = read.delim("data/data1.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]

pA = ggplot(data1,aes(y=expressed_overused_background_WC_duet_ambiguous,x=clade_group,fill=clade_group,label=species))  +
  geom_hline(size=1,linetype="dashed",col="red", yintercept = 0 ) +
  geom_boxplot(alpha=.1) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica",margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0,angle = 50, family="economica"),
    title =  element_text(color="black", size=0, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +ylab("Difference in proportion of WCp codons between\nthe top 5% and bottom 50% expressed")+
  xlab("")  + scale_y_continuous(breaks = seq(-10,25,5))

pA = ggMarginal(pA, type="histogram",fill=set_color[1])
pA

jpeg(paste(path_pannel,"F6pA.jpg",sep=""), width = 5500/1, height = 3000/1,res=460/1)
print(pA)
dev.off()


# Pannel 6 B
dt_graph = data1
ylabel = "expressed_overused_background_WC_duet_ambiguous"
xlabel = "expressed_overused_background_WB_WC_notambiguous"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)

pB =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species"))  +
  geom_abline() +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=26,family="economica"),
    axis.title.y = element_text(color="black", size=26,hjust=c(1,1), family="economica",margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.y =  element_text(color="black", size=24, family="economica"),
    axis.text.x =  element_text(color="black", size=24, family="economica"),
    title =  element_text(color="black", size=20, family="economica"),
    text =  element_text(color="black", size=31, family="economica"),
    legend.text =  element_text(color="black", size=24, family="economica",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.7, face= "italic", size=20, family="economica"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste("LM: "," R"^2,lm_eqn," / PGLS:"," R"^2,pgls_eq), list(nbspecies=nrow(dt_graph),
                                                                                     lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
                                                                                     pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
    # title = substitute(paste("Nspecies = ",nbspecies), list(nbspecies=nrow(dt_graph),
    #                                                         lm_eqn=lm_eqn(lm(lm_y ~ lm_x)),
    #                                                         pgls_eq=lm_eqn(pgls(pgls_y~pgls_x,shorebird))))
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  xlab("Difference in proportion of optimal codons between\nthe top 5% and bottom 50% expressed (%)")  + 
  ylab("Difference in proportion of WCp codons between\nthe top 5% and bottom 50% expressed        ") + scale_y_continuous(breaks = seq(-10,25,5))

pB

jpeg(paste(path_pannel,"F6pB.jpg",sep=""),width = 5200/2, height = 4000/2,res=600/2)
print(pB)
dev.off()


# Figure 6

imgA = load.image(paste(path_pannel,"F6pA.jpg",sep="") )
imgB = load.image(paste(path_pannel,"F6pB.jpg",sep="") )
clade_png<-readPNG(paste(path_require,"clade.png",sep=""))

{
  pdf(file= paste(path_figure,"Figure6.pdf",sep=""), width=7, height=7)
  
  m = matrix(rep(NA,100*9), nrow=9)
  
  for(i in 1:5){
    m[i,]=c(rep(1,10))
  }
  for(i in 5:9){
    m[i,]=c(rep(2,10))
  }
  layout(m)
  
  par(mar=c(1,0, 1, 7))
  plot(imgA, axes=FALSE)
  mtext("A",at=200,adj=-1.5, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(1,0, 0, 0))
  xmonkey=5500
  ymonkey=500
  rasterImage(clade_png,xleft=0+xmonkey, ybottom=800/.4+ymonkey, xright=500/.4+xmonkey, ytop=ymonkey)
  
  par(mar=c(3,0, 2, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=0,adj=-5, side=2, line=1, font=2, cex=1.7,las=2)
  dev.off()
}

