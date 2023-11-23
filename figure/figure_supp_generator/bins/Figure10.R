source("figure/figure_main generator/library_path.R")


############## Pannel 9 A

Freq_opti = read.delim("data/Freq_opti.tab")
Freq_opti[,c("clade_group","lifespan","length","weight")] = clade_dt[Freq_opti$species,c("clade_group","lifespan","length","weight")]

arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
list_species = arbrePhylo$tip.label
Freq_opti = Freq_opti[Freq_opti$species %in% list_species,]

dt_graph = Freq_opti[Freq_opti$type_aa == "WC_duet_ambiguous",]
dt_graph = dt_graph[dt_graph$species %in% clade_dt$species & dt_graph$clade_group %in% c("Lepido Diptera"  ,      "Hymenoptera"    ,   "Other Insecta" ),]
# dt_graph = dt_graph[dt_graph$species %in% clade_dt$species ,]
dt_graph = dt_graph[ dt_graph$nb_codon_not_decoded == 0 & dt_graph$nb_genes > 5000 & dt_graph$pval_aa_fpkm < 0.05  ,]
dt_graph$ecart = (dt_graph$optifreq_top5-dt_graph$opti_freq_low50) - (dt_graph$optifreq_top5_intron-dt_graph$opti_freq_low50_intron)
dt_graph$ecart = dt_graph$optifreq_top5-dt_graph$opti_freq_low50

ylabel = "optifreq_top5"
xlabel = "optifreq_top5_intron"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]),] 
lm_y = dt_graph[,ylabel]
lm_x = dt_graph[,xlabel]
m1 <- lm(lm_y~lm_x)  #Create a linear model

new_x_label = "optifreq_top5_intron"
new_x = dt_graph[,new_x_label]
predicted_y <- predict(m1, newdata = data.frame(x = new_x))
new_y_label = "optifreq_top5"
new_y = dt_graph[,new_y_label]
dt_graph$residual <- new_y - predicted_y


shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)

gls = GLS(shorebird)

p9A = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel),aes(fill=clade_group)) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=20,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black", size=20, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=15, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  # ggtitle(paste("Wb_WC_notambiguous N=",nrow(dt_graph)," / LM:",lm_eqn(lm(lm_y ~ lm_x)),
  #               "\n            PGLS:",lm_eqn(pgls(pgls_y~pgls_x,shorebird)),
  #               "\n  Best",gls[[3]],":",lm_eqn(gls[[2]])))+ annotation_logticks(sides="b") +
  geom_smooth(data = dt_graph,aes(fill="tout"), method = lm, se = FALSE)
p9A




get_CM_dNdS_true<-function(D) {
  # Compute the cumulated number of substitutions over all genes
  cum_KS=sum(D$num_dS)
  cum_KN=sum(D$num_dN)
  cum_OS=sum(D$den_dS/D$branch_length)
  cum_ON=sum(D$den_dN/D$branch_length)
  
  # Compute cumulated DN, DS
  cum_dS=cum_KS/cum_OS
  cum_dN=cum_KN/cum_ON
  
  # Compute cumulated dN/dS
  cum_dNdS=cum_dN/cum_dS
  
  return(cum_dNdS)
}

data_dNdS_euk_v7 = read.delim(paste("/home/fbenitiere/data/Projet-SplicedVariants/DnDs/Metazoa_v9/subset_200_ksites_GC3/data_calculation.tab",sep=""))
data_dNdS_euk_v7 = data_dNdS_euk_v7[!is.na(data_dNdS_euk_v7$branch_length),]

dt_graph$dnds = NA
for (species in unique(dt_graph$species)){
  data_sequence = data_dNdS_euk_v7[data_dNdS_euk_v7$species == species,]
  dt_graph[dt_graph$species == species,]$dnds = get_CM_dNdS_true(data_sequence)
}


ylabel = "residual"
xlabel = "dnds"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]),]
lm_y = dt_graph[,ylabel]
lm_x = (dt_graph[,xlabel])
shorebird <- comparative.data(arbrePhylo, 
                              data.frame(species=dt_graph$species,
                                         pgls_x=lm_x,
                                         pgls_y=lm_y), species, vcv=TRUE)

gls = GLS(shorebird)


p9A = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel),aes(fill=clade_group)) + 
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=20,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",size=20, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=15, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) +
  # scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,365,3650,36500),labels=c(0.05,0.1,0.5,1,5,10,100,"1 yrs","10 yrs","100 yrs")) + xlab("Longevity (days log scale)")+ annotation_logticks(sides="b") +
  ggtitle(paste("N=",nrow(dt_graph)," / LM:",lm_eqn(lm(lm_y ~ lm_x)),
                "\n            PGLS:",lm_eqn(pgls(pgls_y~pgls_x,shorebird)),
                "\n  Best",gls[[3]],":",lm_eqn(gls[[2]])))
p9A

