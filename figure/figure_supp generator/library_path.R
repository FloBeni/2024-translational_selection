
options( stringsAsFactors = F, scipen = 999 )
library(seqinr)
library(stringr)
library(ape)
library(png)
library(ggtree)
library(caper)
library(ggExtra)
library(phylolm)
library(imager)
library(ggplot2)
library(RColorBrewer)
library(stringi)
library(forcats)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))

path_figure = "figure/figure_supp/"
path_pannel = "figure/pannels/"
path_require = "figure/images_library/"


wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")

Clade_color = c("Other Invertebrates"="#f5b48a","Lepido Diptera"="red","Other Tetrapodes"="#A6CEE3","Other Insecta"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Hymenoptera="#ba8e18",Aves="#5b5b5b",Mammalia="#66281A",Embryophyta="#33A02C"
)

Clade_color = Clade_color[c("Embryophyta","Lepido Diptera","Hymenoptera",
                            "Other Insecta","Nematoda","Other Invertebrates","Teleostei",
                            "Mammalia","Aves","Other Tetrapodes")]

arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))


clade_dt = read.delim(paste( "data/clade_dt.tab",sep=""),header=T)
rownames(clade_dt) = clade_dt$species
clade_dt$clade_group = factor(clade_dt$clade_group, levels = c("Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Teleostei","Mammalia","Aves","Other Tetrapodes"))


listNomSpecies = tapply(clade_dt$species,clade_dt$clade_group,function(x)  str_replace_all(x,"_"," "))

# clade_dt = clade_dt[clade_dt$species%in%list_species,]
table(clade_dt$clade_group)

lm_eqn <- function(m=lm(Y ~ X,data)){
  # paste("R2 =", format(summary(m)$r.squared, digits = 2) , "; p-value = ",formatC(summary(m)$coefficients[2,4], format = "e", digits = 2))
  paste("R2 =", round(summary(m)$r.squared, 2) , "; p-value = ",formatC(summary(m)$coefficients[2,4], format = "e", digits = 0))
}


GLS <- function(dataframe=shorebird){
  aic = 1000000
  dt = data.frame()
  for (model in c("LM","lambda","OUfixedRoot","OUrandomRoot","BM")){
    for (measurement_error in c(T,F)){
      if (model == "LM"){
        fit = lm(pgls_y~pgls_x, data = dataframe$data)
        measurement_error = NA
      } else if (model != "lambda"){
        fit <- phylolm(pgls_y~pgls_x, phy = dataframe$phy, data = dataframe$data, model = model,measurement_error=measurement_error)
      } else{ fit <- phylolm(pgls_y~pgls_x, phy = dataframe$phy, data = dataframe$data, model = model)
      measurement_error = NA}
      a = summary(fit)
      if (length(a$optpar)==0){a$optpar=NA}
      if (length(a$aic)==0){a$aic=NA
      a$logLik=NA
      a$optpar=NA
      a$sigma2=NA}
      
      dt = rbind(dt,data.frame(
        model,
        measurement_error,
        p_val_slope = a$coefficients[2,4],
        r.squared = a$r.squared,
        adj.r.squared = a$adj.r.squared,
        aic = a$aic,
        logLik = a$logLik,
        optpar = a$optpar,
        sigma2 = a$sigma2
      ))
      if ( !is.na(a$aic < aic) & a$aic < aic ){ best_fit_model = fit
      best_model = model
      aic = a$aic}
    }
  }
  dt = dt[!duplicated(dt$aic),]
  dt = dt[order(dt$aic),]
  return(list(dt,best_fit_model,best_model))
}
