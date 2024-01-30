# Generate Data 12
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
data1 = data1[data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000,]

tRNA_abundance = tRNA_abundance[rownames(tRNA_abundance) %in% data1$species,]

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon

code = code[!code$aa_name %in% c("Ter"),]

rownames(code) = code$anticodon

data12 = data.frame()
for (anticodon in code$anticodon){
  dt = data.frame(abundance = tRNA_abundance[,anticodon])
  dt$species = rownames(tRNA_abundance)
  dt$nb_syn = code[anticodon,]$nb_syn
  dt$amino_acid = code[anticodon,]$aa
  dt$anticodon = anticodon
  dt$codon = code[anticodon,]$codon
  data12 = rbind(data12,dt)
}

for (species in GTDrift_list_species$species){
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  rownames(tRNA_optimal) = tRNA_optimal$codon
  data12[data12$species == species,c("POC1","POC2")] = tRNA_optimal[ data12[data12$species == species,]$codon,c("POC1","POC2")]
}  



data12$color = sapply(data12$codon,function(x)substr(x,3,3))

data12$amino_acid = factor(data12$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa))
data12$anticodon = str_replace_all(data12$anticodon,'T','U')
data12$codon = str_replace_all(data12$codon,'T','U')

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
data12$title = paste(data12$anticodon,"(",data12$codon,")",sep="")
data12$codon = factor(data12$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
data12$title = factor(data12$title,levels= tapply(data12$title, as.integer(data12$codon),unique))

nb_sp = length(unique(data12$species))
data12$nb_species_0 = tapply(data12$abundance != 0,data12$codon,sum)[data12$codon]
data12$nb_species_0 = round(data12$nb_species_0 / nb_sp*100)
data12$y_axis_0 = tapply(data12$abundance ,data12$codon,function(x) quantile(x,0.9))[data12$codon]
data12[duplicated(data12$codon) | data12$nb_species_0 > 30,]$nb_species_0 = NA
data12[!is.na(data12$nb_species_0),]$nb_species_0 = paste(data12[!is.na(data12$nb_species_0),]$nb_species_0 ,"%")


write.table(data12,"data/data12_supp.tab",quote=F,row.names = F,sep="\t")


