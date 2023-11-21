
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
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))

path_figure = "figure/figure_main/"
path_pannel = "figure/pannels/"
path_require = "figure/images_library/"


wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")

Clade_color = c("Other Invertebrates"="#f5b48a","Mecopterida"="red","Other Tetrapods"="#A6CEE3","Other Insecta"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Hymenoptera="#ba8e18",Aves="#5b5b5b",Mammalia="#66281A",Embryophyta="#33A02C"
)

Clade_color = Clade_color[c("Embryophyta","Mecopterida","Hymenoptera",
                            "Other Insecta","Nematoda","Other Invertebrates",
                            "Mammalia","Aves","Other Tetrapods","Teleostei")]

arbrePhylo = read.tree(paste("data/GTDrift_metazoa_phylogenetic_tree.nwk",sep=""))

life_history_traits = read.delim("data/GTDrift_life_history_traits.tab")
rownames(life_history_traits) = paste(life_history_traits$species,life_history_traits$life_history_traits,sep="_")
  
GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species
GTDrift_list_species[GTDrift_list_species$clade_group == "Other Vertebrates" ,]$clade_group = "Other Tetrapods"

GTDrift_list_species$clade_group = factor(GTDrift_list_species$clade_group, levels = c("Mecopterida","Hymenoptera","Other Insecta",
                                                                                       "Nematoda","Other Invertebrates","Teleostei",
                                                                                       "Mammalia","Aves","Other Tetrapods"))

GTDrift_list_species$length_cm = life_history_traits[paste(GTDrift_list_species$species,"length_cm",sep="_"),]$value
GTDrift_list_species$lifespan_days = life_history_traits[paste(GTDrift_list_species$species,"lifespan_days",sep="_"),]$value
GTDrift_list_species$weight_kg = life_history_traits[paste(GTDrift_list_species$species,"weight_kg",sep="_"),]$value


listNomSpecies = tapply(GTDrift_list_species$species,GTDrift_list_species$clade_group,function(x)  str_replace_all(x,"_"," "))



lm_eqn <- function(m=lm(Y ~ X,data)){
  pvalue = summary(m)$coefficients[2,4]
  if (pvalue<10^-100){
    paste(" = ", round(summary(m)$r.squared, 2) , ", p-value = 0",sep="")
  } else {
    paste(" = ", round(summary(m)$r.squared, 2) , ", p-value = ",formatC(pvalue, format = "e", digits = 0),sep="")
  }
}




