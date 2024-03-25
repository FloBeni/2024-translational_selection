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
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
data1 = data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species != "Eumeta_japonica" & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 & data1$nb_codon_not_decoded == 0,]

tRNA_abundance = tRNA_abundance[rownames(tRNA_abundance) %in% data1$species,]

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon

code = code[!code$aa_name %in% c("Ter"),]

rownames(code) = code$anticodon

data14 = data.frame()
for (anticodon in code$anticodon){
  dt = data.frame(abundance = tRNA_abundance[,anticodon])
  dt$species = rownames(tRNA_abundance)
  dt$nb_syn = code[anticodon,]$nb_syn
  dt$amino_acid = code[anticodon,]$aa
  dt$anticodon = anticodon
  dt$codon = code[anticodon,]$codon
  data14 = rbind(data14,dt)
}

for (species in GTDrift_list_species$species){
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  rownames(tRNA_optimal) = tRNA_optimal$codon
  data14[data14$species == species,c("POC1","POC2")] = tRNA_optimal[ data14[data14$species == species,]$codon,c("POC1","POC2")]
}  



data14$color = sapply(data14$codon,function(x)substr(x,3,3))

data14$amino_acid = factor(data14$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa))
data14$anticodon = str_replace_all(data14$anticodon,'T','U')
data14$codon = str_replace_all(data14$codon,'T','U')

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
data14$title = paste(data14$anticodon,"(",data14$codon,")",sep="")
data14$codon = factor(data14$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
data14$title = factor(data14$title,levels= tapply(data14$title, as.integer(data14$codon),unique))

nb_sp = length(unique(data14$species))
data14$nb_species_0 = tapply(data14$abundance == 0,data14$codon,sum)[data14$codon]
data14$nb_species_0 = round(data14$nb_species_0 / nb_sp*100)
data14$y_axis_0 = tapply(data14$abundance ,data14$codon,function(x) quantile(x,0.9))[data14$codon]
data14[duplicated(data14$codon) ,]$nb_species_0 = NA
data14[!is.na(data14$nb_species_0),]$nb_species_0 = paste(data14[!is.na(data14$nb_species_0),]$nb_species_0 ,"%")


write.table(data14,"data/data14_supp.tab",quote=F,row.names = F,sep="\t")


