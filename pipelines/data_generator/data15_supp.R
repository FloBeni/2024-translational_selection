# Generate Data 15

options(scipen=999)

library(stringr)
library(ggplot2)


path = "/home/fbenitiere/data/"
# path = "/beegfs/data/fbenitiere/"

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab",comment.char = "#")
rownames(GTDrift_list_species) = GTDrift_list_species$species

data15 = data.frame()
for (species in GTDrift_list_species$species){print(species)
    dt = read.delim(paste(path,"Projet-NeGA/translational_selection/GC_gap_per_window_per_species/",species,".tab",sep=""),comment.char = "#")
    
    size_windows = 100
    value_windows = floor(seq(0,max(dt$start),100)/size_windows) *size_windows+size_windows/2
    names(value_windows) = as.character(seq(0,max(dt$start),100))


    dt$group = value_windows[as.character(dt$start)]

    dg = data.frame(
      species,
      per_windows=names(tapply(dt$pos3_sites,paste(dt$from,dt$group),sum)),
      nb_genes=tapply(dt$protein,paste(dt$from,dt$group),function(x) length(unique(x))),
      pos3_sites = tapply(dt$pos3_sites,paste(dt$from,dt$group),sum),
      gc3_count = tapply(dt$gc3_count,paste(dt$from,dt$group),sum),
      posi_sites = tapply(dt$posi_sites,paste(dt$from,dt$group),sum),
      gci_count = tapply(dt$gci_count,paste(dt$from,dt$group),sum)
    )

    data15 = rbind(data15,dg)

}

write.table(data15,"data/data15_supp.tab",quote=F,row.names = F,sep="\t")

