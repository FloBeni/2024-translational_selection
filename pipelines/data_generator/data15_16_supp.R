

library(stringi)

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code = code[code$aa_name != "Ter",]



tRNA_abundance = read.delim("/home/fbenitiere/2024-translational_selection/data/tRNA_abundance.tab")

species="Blattella_germanica"
dt = data.frame(nb_copy=tapply(code$anticodon,code$aa_name,function(x) sum(tRNA_abundance[species,x])))
code$nb_copy = unlist(tRNA_abundance[species,code$anticodon])
code[code$anticodon == "AGA",]
min(code[code$anticodon != "AGA" & code$nb_copy!=0,]$nb_copy)
max(code[code$anticodon != "AGA" & code$nb_copy!=0,]$nb_copy)

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

all_data = data.frame()
for (species in GTDrift_list_species$species){
  
  dt = data.frame(nb_copy=tapply(code$anticodon,code$aa_name,function(x) sum(tRNA_abundance[species,x])))
  dt$aa_name = rownames(dt)
  dt$species = species
  
  all_data = rbind(all_data,dt)
  
}

data16 = cbind(all_data[all_data$species == "Caenorhabditis_elegans", ],all_data[all_data$species == "Hydra_vulgaris", ])
colnames(data16) = c("nb_copy_Caenorhabditis_elegans" ,"aa_name_Caenorhabditis_elegans", "species", "nb_copy_Hydra_vulgaris" ,"aa_name_Hydra_vulgaris" ,"species")

write.table(data16,"data/data_16.tab",quote=F,row.names = F,sep="\t")

all_cor = data.frame()
species1 = "Hydra_vulgaris"
for (species2 in GTDrift_list_species$species){
  
  spearman_method_aa = cor.test( all_data[all_data$species == species1, ]$nb_copy, all_data[all_data$species == species2, ]$nb_copy,method="spearman",exact=F)
  
  all_cor = rbind(all_cor,data.frame(
    species1,
    species2,
    rho = spearman_method_aa$estimate,
    pval = spearman_method_aa$p.value
  ))
  
}


write.table(all_cor,"data/data_15.tab",quote=F,row.names = F,sep="\t")

