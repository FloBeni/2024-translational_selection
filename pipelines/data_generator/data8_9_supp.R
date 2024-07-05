# Generate Data 8 and 9
library(stringi)

code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
rownames(code) = code$codon
code = code[code$aa_name != "Ter",]

data7 <- read.table("data/data7_supp.tab",header=T)

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab",comment.char = "#")
rownames(GTDrift_list_species) = GTDrift_list_species$species

all_data = data.frame()
for (species in GTDrift_list_species$species){
  
  dt = data.frame(nb_copy=tapply(code$anticodon,code$aa_name,function(x) sum(data7[species,x])))
  dt$aa_name = rownames(dt)
  dt$species = species
  
  all_data = rbind(all_data,dt)
  
}

data9 = cbind(all_data[all_data$species == "Caenorhabditis_elegans", ],all_data[all_data$species == "Hydra_vulgaris", ])
colnames(data9) = c("nb_copy_Caenorhabditis_elegans" ,"aa_name_Caenorhabditis_elegans", "species", "nb_copy_Hydra_vulgaris" ,"aa_name" ,"species")
data9 = data9[,c("aa_name","nb_copy_Caenorhabditis_elegans" , "nb_copy_Hydra_vulgaris" )]
write.table(data9,"data/data9_supp.tab",quote=F,row.names = F,sep="\t")

data8 = data.frame()
species1 = "Hydra_vulgaris"
for (species2 in GTDrift_list_species$species){
  
  spearman_method_aa = cor.test( all_data[all_data$species == species1, ]$nb_copy, all_data[all_data$species == species2, ]$nb_copy,method="spearman",exact=F)
  
  data8 = rbind(data8,data.frame(
    species1,
    species2,
    rho = spearman_method_aa$estimate,
    pval = spearman_method_aa$p.value
  ))
  
}


write.table(data8,"data/data8_supp.tab",quote=F,row.names = F,sep="\t")

