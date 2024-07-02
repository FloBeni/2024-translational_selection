# Generate Data 4
options(stringsAsFactors = F, scipen = 999)
library(stringr)
library(stringi)

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

### tRNA abundance
tRNA_abundance = data.frame()
for (species in GTDrift_list_species$species){
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  dt = t(tRNA_optimal[,c("anticodon","nb_tRNA_copies")])
  dt = data.frame(dt)
  colnames(dt) = dt[1,]  
  dt = dt[2,]
  rownames(dt) = species
  tRNA_abundance = rbind(tRNA_abundance,dt)
}  
tRNA_abundance <- data.frame(sapply( tRNA_abundance, as.numeric ))
rownames(tRNA_abundance) = GTDrift_list_species$species

data1 = read.delim("data/data1_supp.tab")
data1 = data1[data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 & data1$nb_codon_not_decoded == 0,]

tRNA_abundance = tRNA_abundance[rownames(tRNA_abundance) %in% data1$species,]

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon

code = code[!code$aa_name %in% c("Ter"),]

rownames(code) = code$anticodon

data4 = data.frame()
for (anticodon in code$anticodon){
  dt = data.frame(abundance = tRNA_abundance[,anticodon])
  dt$species = rownames(tRNA_abundance)
  dt$nb_syn = code[anticodon,]$nb_syn
  dt$amino_acid = code[anticodon,]$aa
  dt$anticodon = anticodon
  dt$codon = code[anticodon,]$codon
  data4 = rbind(data4,dt)
}

for (species in GTDrift_list_species$species){
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  rownames(tRNA_optimal) = tRNA_optimal$codon
  data4[data4$species == species,c("POC1","POC2")] = tRNA_optimal[ data4[data4$species == species,]$codon,c("POC1","POC2")]
}  



data4$color = sapply(data4$codon,function(x)substr(x,3,3))

data4$amino_acid = factor(data4$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa))
data4$anticodon = str_replace_all(data4$anticodon,'T','U')
data4$codon = str_replace_all(data4$codon,'T','U')

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
data4$title = paste(data4$anticodon,"(",data4$codon,")",sep="")
data4$codon = factor(data4$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
data4$title = factor(data4$title,levels= tapply(data4$title, as.integer(data4$codon),unique))

nb_sp = length(unique(data4$species))
data4$nb_species_0 = tapply(data4$abundance == 0,data4$codon,sum)[data4$codon]
data4$nb_species_0 = round(data4$nb_species_0 / nb_sp*100)
data4$y_axis_0 = tapply(data4$abundance ,data4$codon,function(x) quantile(x,0.9))[data4$codon]
# data4[duplicated(data4$codon) | data4$nb_species_0 < 50,]$nb_species_0 = NA
data4[duplicated(data4$codon) ,]$nb_species_0 = NA
data4[!is.na(data4$nb_species_0),]$nb_species_0 = paste(data4[!is.na(data4$nb_species_0),]$nb_species_0 ,"%")


write.table(data4,"data/data4_supp.tab",quote=F,row.names = F,sep="\t")


