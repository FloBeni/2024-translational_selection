source("figure/figure_main generator/library_path.R")


############## Pannel 7 A

data12 = read.delim("data/data12.tab")
data12$clade_group = clade_dt[data12$species,]$clade_group

data12 = data12[ data12$nb_codon_not_decoded == 0 & data12$pval_aa_fpkm < 0.05 ,]

data12$ecart = (data12$optifreq_top5-data12$opti_freq_low50) - (data12$optifreq_top5_intron-data12$opti_freq_low50_intron)
data12$ecart = data12$ecart * 100

data12[,c("clade_group","lifespan","length","weight")] = clade_dt[data12$species,c("clade_group","lifespan","length","weight")]


arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label
data12 = data12[data12$species %in% list_species,]

dnds = read.delim(paste("/home/fbenitiere/2024-EukGTDrift/data/dnds_0.1_1_dS.tab",sep=""))
rownames(dnds) = dnds$species
dnds[dnds$species %in% c("Cervus_elaphus","Vulpes_lagopus","Bombus_vancouverensis","Gymnodraco_acuticeps"),]$dNdS = NA
data12$dnds = dnds[data12$species,]$dNdS

dt_graph = data12[data12$type_aa == "Wb_WC_notambiguous",]
colnames(dt_graph)
ylabel = "ecart"
xlabel = "var_gci_top5"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = (dt_graph[,xlabel])
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)


p7A = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) + ggtitle("Wb_WC_notambiguous") +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=.8) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=16, family="economica"),
    legend.text =  element_text(color="black", size=15, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("Difference in proportion of OC between\nthe top 5% and bottom 50% expressed (%)") + 
  # scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000),labels=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000)) + xlab("Longevity (days log scale)")+
  xlab("Variance per gene GCi\namong top 5% expressed genes")+
  # xlab("dN/dS")+
  ggtitle(paste(
    # "Wb_WC_notambiguous N=",nrow(dt_graph),
    "LM: ",lm_eqn(lm(lm_y ~ lm_x)),
    " / PGLS: ",lm_eqn(pgls(pgls_y~pgls_x,shorebird))
    ,sep="")) + ylim(-5,20)
+ annotation_logticks(sides="b") 
p7A
ggplotly(p7A)
jpeg(paste(path_pannel,"p7A.jpg",sep=""), width = 4000/1, height = 3000/1,res=550/1)
print(p7A)
dev.off()



############## Pannel 7 B
dt_graph = data12[data12$type_aa == "Wb_WC_notambiguous",]

ylabel = "ecart"
xlabel = "length"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
lm_y = dt_graph[,ylabel]
lm_x = log10(dt_graph[,xlabel])
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)

p7B = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group")) + ggtitle("WC_duet_ambiguous") +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=15, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  ylab("") + 
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000,5000),labels=c(0.01,0.1,1,10,100,1000,5000)) + xlab("Body length (cm log scale)")+
  ggtitle(paste(
    "LM: ",lm_eqn(lm(lm_y ~ lm_x)),
    " / PGLS: ",lm_eqn(pgls(pgls_y~pgls_x,shorebird)),sep=""
  ))+ annotation_logticks(sides="b")+ ylim(-5,20)
p7B

jpeg(paste(path_pannel,"p7B.jpg",sep=""), width = 4000/1, height = 3000/1,res=600/1)
print(p7B)
dev.off()




############## Pannel 7 C
data13 = read.delim("data/data13.tab")
rownames(data13) = data13$species

data12$GC_DUC = data13[data12$species,]$GC_DUC
data12$GC_DUC_quart = data13[data12$species,]$GC_DUC_quart
data12$prop_abond_quart = data13[data12$species,]$prop_abond_quart
data12$prop_abond = data13[data12$species,]$prop_abond
data12$mean_ecart = data13[data12$species,]$mean_ecart
data12$mean_ecart_quart = data13[data12$species,]$mean_ecart_quart

dt_graph = data12[data12$type_aa == "Wb_WC_notambiguous",]
dt_graph = dt_graph[dt_graph$clade_group %in% c( "Lepido Diptera"),]
# dt_graph = dt_graph[dt_graph$species != "Eumeta_japonica",]
# dt_graph = dt_graph[dt_graph$ecart > 5 , ]

{
  ylabel = "ecart"
  xlabel = "prop_abond_quart"
  lm_y = dt_graph[,ylabel]
  lm_x = dt_graph[,xlabel]
  shorebird <- comparative.data(arbrePhylo, 
                                data.frame(species=dt_graph$species,
                                           pgls_x=lm_x,
                                           pgls_y=lm_y), species, vcv=TRUE)
  
  gls = GLS(shorebird)
  dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]),] 
  
  p7C = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
    # geom_errorbar(aes(ymin = gc3-std_gc3, ymax = std_gc3+gc3),alpha=0.2) +
    # geom_errorbarh(aes(xmin = gci-std_gci, xmax = std_gci+gci),alpha=0.2) +
    geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
      axis.title.x = element_text(color="black", size=25,family="economica"),
      axis.title.y = element_text(color="black", size=25, family="economica"),
      axis.text.y =  element_text(color="black", size=20, family="economica"),
      axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
      title =  element_text(color="black", size=15, family="economica"),
      legend.text =  element_text(color="black", size=15, family="economica")
    ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
    ylab("Difference in proportion of OC between\nthe top 5% and bottom 50% expressed") + 
    xlab("Proportion of PC that are OC") +
    ggtitle(paste("N=",nrow(dt_graph),"\nLM:",lm_eqn(lm(lm_y ~ lm_x))
                  # "\nPGLS:",lm_eqn(pgls(pgls_y~pgls_x,shorebird))
                  # , "\n  Best",gls[[3]],":",lm_eqn(gls[[2]])
    ))
  print(p7C)
}
# library(plotly)
# ggplotly(p7C)

jpeg(paste(path_pannel,"p7C.jpg",sep=""), width = 4000/1, height = 3000/1,res=500/1)
print(p7C)
dev.off()


############## Pannel 7 D
colnames(dt_graph)

# dt_graph = dt_graph[dt_graph$ecart > 5 , ]
{
  ylabel = "GC_DUC_quart"
  xlabel = "gci"
  lm_y = dt_graph[,ylabel]
  lm_x = dt_graph[,xlabel]
  shorebird <- comparative.data(arbrePhylo, 
                                data.frame(species=dt_graph$species,
                                           pgls_x=lm_x,
                                           pgls_y=lm_y), species, vcv=TRUE)
  
  # gls = GLS(shorebird)
  dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]),] 
  
  # dt_graph[,xlabel] = rank(dt_graph[,c(xlabel)])
  # dt_graph[,ylabel] = rank(dt_graph[,c(ylabel)])
  
  p7D = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
    # geom_errorbar(aes(ymin = gc3-std_gc3, ymax = std_gc3+gc3),alpha=0.2) +
    # geom_errorbarh(aes(xmin = gci-std_gci, xmax = std_gci+gci),alpha=0.2) +
    geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
      axis.title.x = element_text(color="black", size=25,family="economica"),
      axis.title.y = element_text(color="black", size=25, family="economica"),
      axis.text.y =  element_text(color="black", size=20, family="economica"),
      axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
      title =  element_text(color="black", size=15, family="economica"),
      legend.text =  element_text(color="black", size=15, family="economica")
    ) + scale_fill_manual(values=Clade_color) +
    ylab("GC3 content of PC") +
    # ylab("Average per gene GC3") +
    xlab("Average per gene GCi") +
    ggtitle(paste("N=",nrow(dt_graph),"\nLM:",lm_eqn(lm(lm_y ~ lm_x))
                  # "\nPGLS:",lm_eqn(pgls(pgls_y~pgls_x,shorebird))
                  # , "\n  Best",gls[[3]],":",lm_eqn(gls[[2]])
    ))+
    guides(fill = FALSE, linetype=F, pch=F, alpha=F)
  print(p7D)
}

# library(plotly)
# ggplotly(p7D)

jpeg(paste(path_pannel,"p7D.jpg",sep=""), width = 4000/1, height = 3000/1,res=500/1)
print(p7D)
dev.off()

############## Pannel 7 D
# data7 = read.delim("data/data7.tab")
# data7$clade_group = clade_dt[data7$species,]$clade_group
# 
# p7D = ggplot(data7[ data7$genome_character == "CG",],aes(y=shift*100,x=clade_group,fill=clade_group)) +
#   geom_hline(yintercept=0,alpha=0.2)+ geom_boxplot() + 
#   scale_fill_manual("Clades",values=Clade_color) + scale_size(range=c(0,8),guide="legend") + theme_bw() + theme(
#     axis.title.x = element_text(color="black", size=25,family="economica"),
#     axis.title.y = element_text(color="black", size=25, family="economica"),
#     axis.text.y =  element_text(color="black", size=20, family="economica"),
#     axis.text.x =  element_text(color="black",vjust=.5, size=0,angle=60, family="economica"),
#     title =  element_text(color="black", size=15, family="economica"),
#     legend.text =  element_text(color="black", size=20, family="economica")
#   ) + ylab("CG Shift observed - expected (%)") + xlab("") + theme(legend.position='none')
# p7D
# 
# jpeg(paste(path_pannel,"p7D.jpg",sep=""), width = 4000/1, height = 3000/1,res=500/1)
# print(p7D)
# dev.off()

############## Figure 7

imgA = load.image(paste(path_pannel,"p7A.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p7B.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p7C.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p7D.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure8.pdf",sep=""), width=7, height=6)
  
  m = matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:5){
    m[,i]=c(rep(1,5),rep(2,5))
  }
  
  for(i in 6:10){
    m[,i]=c(rep(3,5),rep(4,5))
  }
  layout(m)
  m
  
  par(mar=c(1, 0, 1.7, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 1.7, 0))
  plot(imgC, axes=FALSE)
  mtext("C",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 1.7, 0))
  plot(imgB, axes=FALSE)
  mtext("B",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(1, 0, 1.7, 1))
  plot(imgD, axes=FALSE)
  mtext("D",at=-50,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
  dev.off()
}

