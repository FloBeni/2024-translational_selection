
options( stringsAsFactors = F, scipen = 999 )
library(seqinr)
library(stringr)
library(ape)
library(patchwork)
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
library(scales)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))

path_figure = "figure/figure_main/"
path_pannel = "figure/pannels/"
path_require = "figure/images_library/"


wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")

Clade_color = c("Other Metazoans"="#f5b48a","Diptera"="red","Other Tetrapods"="#A6CEE3","Other Insecta"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Hymenoptera="#ba8e18",Aves="#5b5b5b",Mammalia="#66281A",Lepidoptera="#33A02C",Coleoptera="#f1dd41"
)

Clade_color = Clade_color[c("Embryophyta","Diptera","Lepidoptera","Coleoptera","Hymenoptera",
                            "Other Insecta","Nematoda","Other Metazoans",
                            "Mammalia","Aves","Other Tetrapods","Teleostei")]



arbrePhylo = read.tree(paste("data/GTDrift_metazoa_phylogenetic_tree.nwk",sep=""))

life_history_traits = read.delim("data/GTDrift_life_history_traits.tab")
rownames(life_history_traits) = paste(life_history_traits$species,life_history_traits$life_history_traits,sep="_")
  
GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

GTDrift_list_species[GTDrift_list_species$clade == "Coleoptera" ,]$clade_group = "Coleoptera"
GTDrift_list_species[GTDrift_list_species$clade == "Diptera" ,]$clade_group = "Diptera"
GTDrift_list_species[GTDrift_list_species$clade == "Lepidoptera" ,]$clade_group = "Lepidoptera"
GTDrift_list_species[GTDrift_list_species$clade_group == "Other Invertebrates" ,]$clade_group = "Other Metazoans"
GTDrift_list_species[GTDrift_list_species$clade_group == "Other Vertebrates" ,]$clade_group = "Other Tetrapods"

GTDrift_list_species$clade_group = factor(GTDrift_list_species$clade_group, levels = c("Diptera","Lepidoptera","Coleoptera","Hymenoptera",
                                                                                        "Other Insecta","Nematoda",
                                                                                        "Mammalia","Aves","Other Tetrapods","Teleostei","Other Metazoans"))

GTDrift_list_species$length_cm = life_history_traits[paste(GTDrift_list_species$species,"length_cm",sep="_"),]$value
GTDrift_list_species$lifespan_days = life_history_traits[paste(GTDrift_list_species$species,"lifespan_days",sep="_"),]$value
GTDrift_list_species$weight_kg = life_history_traits[paste(GTDrift_list_species$species,"weight_kg",sep="_"),]$value


listNomSpecies = tapply(GTDrift_list_species$species,GTDrift_list_species$clade_group,function(x)  str_replace_all(x,"_"," "))



fitted_model <- function(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree = NA){
  dt_fit = data.frame()
  if ( length(tree) != 1){
    shorebird <- comparative.data(tree, 
                                  data.frame(label=label,
                                             x=x,
                                             y=y), label, vcv=TRUE)
    fit = pgls(y~x,shorebird)
    summ_fit = summary(fit)
    dt_fit = rbind(dt_fit,data.frame(
      model="PGLS",
      p_val_slope = summ_fit$coefficients[2,4],
      r.squared = summ_fit$r.squared,
      adj.r.squared = summ_fit$adj.r.squared,
      aic = AIC(fit),
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    ))
  } 
  
  fit = lm(y~x)
  summ_fit = summary(fit)
  dt_fit = rbind(dt_fit,data.frame(
    model="LM",
    p_val_slope = summ_fit$coefficients[2,4],
    r.squared = summ_fit$r.squared,
    adj.r.squared = summ_fit$adj.r.squared,
    aic = AIC(fit),
    slope = coef(fit)[2],
    intercept = coef(fit)[1]
  ))
  
  model_sub = dt_fit[dt_fit$aic != min(dt_fit$aic),]
  dt_fit = dt_fit[dt_fit$aic == min(dt_fit$aic),]
  
  model = paste(dt_fit$model,sep="")
  AIC = paste(round(dt_fit$aic),sep="")
  R2 = paste(round(dt_fit$r.squared, 2),sep="")
  pvalue = paste(formatC(dt_fit$p_val_slope, format = "e", digits = 0),sep="")
  model_non_opti = ""
  
  if ( length(tree) != 1){
    model_non_opti = paste("/ ",model_sub$model,": AIC = ",round(model_sub$aic),sep="")
  }
  return(list(model=model,aic=AIC,r2=R2,pvalue=pvalue,model_non_opti=model_non_opti,slope=dt_fit$slope,intercept=dt_fit$intercept))
} 


