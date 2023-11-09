

path = "/home/fbenitiere/data/Projet-SplicedVariants/"

clade_dt = rbind( read.table(paste( path,"Fichiers-data/metazoa_species_clade_lht.tab",sep=""),header=T) ,
                  read.table(paste( path,"Fichiers-data/metazoa_species2_clade_lht.tab",sep=""),header=T))
rownames(clade_dt) = clade_dt$species

table(clade_dt$clade)

arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
clade_dt = clade_dt[clade_dt$species %in% arbrePhylo$tip.label,]

df = data.frame()
for (species_name in clade_dt$species){print(species_name)
  table_ncbi = read.delim(paste(path ,"Annotations/" , species_name , "/taxonomy_ncbi.tab" , sep = ""))
  table_ncbi$species = species_name
  df = rbind(df,table_ncbi)
}

clade_dt$clade_group = "Other Invertebrates"
clade_dt[ df[df$name == "Vertebrata",]$species,]$clade_group = "Other Vertebrates"
clade_dt[ df[df$name == "Tetrapoda",]$species,]$clade_group = "Other Tetrapodes"
clade_dt[ df[df$name == "Insecta",]$species,]$clade_group = "Other Insecta"
clade_dt[ df[df$name %in% c("Diptera","Lepidoptera"),]$species,]$clade_group = "Lepido Diptera"
clade_dt[ df[df$name %in% c("Nematoda","Hymenoptera","Mammalia","Aves","Teleostei","Embryophyta"),]$species,]$clade_group =
  clade_dt[ df[df$name %in% c("Nematoda","Hymenoptera","Mammalia","Aves","Teleostei","Embryophyta"),]$species,]$clade


table(clade_dt$clade_group)



write.table(clade_dt,"data/clade_dt.tab",quote=F,row.names = F,sep="\t")
