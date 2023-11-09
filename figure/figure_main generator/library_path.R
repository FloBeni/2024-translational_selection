
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

path_figure = "figure/figure_main/"
path_pannel = "figure/pannels/"
path_require = "figure/images_library/"


wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")

Clade_color = c("Other Invertebrates"="#f5b48a","Lepido Diptera"="red","Other Tetrapods"="#A6CEE3","Other Insecta"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Hymenoptera="#ba8e18",Aves="#5b5b5b",Mammalia="#66281A",Embryophyta="#33A02C"
)

Clade_color = Clade_color[c("Embryophyta","Lepido Diptera","Hymenoptera",
                            "Other Insecta","Nematoda","Other Invertebrates",
                            "Mammalia","Aves","Other Tetrapods","Teleostei")]

arbrePhylo = read.tree(paste("data/GTDrift_metazoa_phylogenetic_tree.nwk",sep=""))



GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species
GTDrift_list_species[GTDrift_list_species$clade_group == "Other Vertebrates" ,]$clade_group = "Other Tetrapods"

GTDrift_list_species$clade_group = factor(GTDrift_list_species$clade_group, levels = c("Lepido Diptera","Hymenoptera","Other Insecta",
                                                                                       "Nematoda","Other Invertebrates","Teleostei",
                                                                                       "Mammalia","Aves","Other Tetrapods"))



listNomSpecies = tapply(GTDrift_list_species$species,GTDrift_list_species$clade_group,function(x)  str_replace_all(x,"_"," "))


lm_eqn <- function(m=lm(Y ~ X,data)){
  paste("R2 = ", round(summary(m)$r.squared, 2) , ", p-value = ",formatC(summary(m)$coefficients[2,4], format = "e", digits = 0),sep="")
}


