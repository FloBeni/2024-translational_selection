# Generate Data 7
library(stringi)

code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
rownames(code) = code$codon

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab",comment.char = "#")
rownames(GTDrift_list_species) = GTDrift_list_species$species


### tRNA abundance
data7 = data.frame()
for (species in GTDrift_list_species$species){
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""),comment.char = "#")
  dt = t(tRNA_optimal[,c("anticodon","nb_tRNA_copies")])
  dt = data.frame(dt)
  colnames(dt) = dt[1,]  
  dt = dt[2,]
  rownames(dt) = species
  data7 = rbind(data7,dt)
}  

write.table(data7,"data/data7_supp.tab",quote=F,row.names = T,sep="\t")