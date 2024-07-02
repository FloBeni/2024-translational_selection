# Get data from GTDrift
options(stringsAsFactors = F, scipen = 999)


# Download data/GTDrift_Metazoa_dNdS.tab and data/GTDrift_Metazoa_phylogenetic_tree.nwk
bash_command <- paste("wget https://zenodo.org/records/10908656/files/dNdS.tar.gz?download=1; tar -xvzf dNdS.tar.gz?download=1; cp database/dNdS/phylogeny/Metazoa_root.nwk data/GTDrift_Metazoa_phylogenetic_tree.nwk; cp database/dNdS/Metazoa.tab data/GTDrift_Metazoa_dNdS.tab; rm -r dNdS.tar.gz?download=1 database",sep="")
system(bash_command)

dnds = read.delim("data/GTDrift_Metazoa_dNdS.tab")


# Download data/GTDrift_list_species.tab and filter
bash_command <- paste("wget https://zenodo.org/records/10908656/files/list_species.tab?download=1; cp list_species.tab?download=1 data/GTDrift_list_species.tab; rm list_species.tab?download=1",sep="")
system(bash_command)
list_species = read.delim("data/GTDrift_list_species.tab")
list_species = list_species[list_species$species %in% dnds$species,]
write.table(list_species,paste("data/GTDrift_list_species.tab",sep=""),sep="\t",quote=F,row.names=F)


# Download data/GTDrift_life_history_traits_and_polymorphism_derived_Ne.tab and filter
bash_command <- paste("wget https://zenodo.org/records/10908656/files/life_history_traits_and_polymorphism_derived_Ne.tab?download=1; cp life_history_traits_and_polymorphism_derived_Ne.tab?download=1 data/GTDrift_life_history_traits_and_polymorphism_derived_Ne.tab; rm life_history_traits_and_polymorphism_derived_Ne.tab?download=1",sep="")
system(bash_command)
list_species = read.delim("data/GTDrift_life_history_traits_and_polymorphism_derived_Ne.tab")
list_species = list_species[list_species$species %in% dnds$species,]
write.table(list_species,paste("data/GTDrift_life_history_traits_and_polymorphism_derived_Ne.tab",sep=""),sep="\t",quote=F,row.names=F)


# Download data/GTDrift_Metazoa_taxonomy.tab
bash_command <- paste("wget https://zenodo.org/records/11034770/files/FloBeni/2024-GTDrift-Revision_NARGAB.zip?download=1; unzip 2024-GTDrift-Revision_NARGAB.zip?download=1 -d GTDrift; cp GTDrift/FloBeni-2024-GTDrift-0f7d2f6/data/taxonomy.tab data/GTDrift_Metazoa_taxonomy.tab; rm -r GTDrift 2024-GTDrift-Revision_NARGAB.zip?download=1",sep="")
system(bash_command)
list_species = read.delim("data/GTDrift_Metazoa_taxonomy.tab")
list_species = list_species[list_species$species %in% dnds$species,]
write.table(list_species,paste("data/GTDrift_Metazoa_taxonomy.tab",sep=""),sep="\t",quote=F,row.names=F)

