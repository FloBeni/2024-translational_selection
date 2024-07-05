# Generate Data 11
options(stringsAsFactors = F, scipen = 999)
library(stringr)
library(stringi)

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab",comment.char = "#")
rownames(GTDrift_list_species) = GTDrift_list_species$species

### tRNA abundance
tRNA_abundance = data.frame()
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
  tRNA_abundance = rbind(tRNA_abundance,dt)
}  
tRNA_abundance <- data.frame(sapply( tRNA_abundance, as.numeric ))
rownames(tRNA_abundance) = GTDrift_list_species$species

data1 = read.delim("data/data1_supp.tab",comment.char = "#")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
data1 = data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species != "Eumeta_japonica" & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 & data1$nb_codon_not_decoded == 0,]

tRNA_abundance = tRNA_abundance[rownames(tRNA_abundance) %in% data1$species,]

code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
rownames(code) = code$codon

code = code[!code$aa_name %in% c("Ter"),]

rownames(code) = code$anticodon

data11 = data.frame()
for (anticodon in code$anticodon){
  dt = data.frame(species = rownames(tRNA_abundance))
  dt$amino_acid = code[anticodon,]$aa
  dt$nb_syn = code[anticodon,]$nb_syn
  dt$anticodon = anticodon
  dt$codon = code[anticodon,]$codon
  dt$tRNA_gene_copy = tRNA_abundance[,anticodon]
  data11 = rbind(data11,dt)
}

for (species in GTDrift_list_species$species){
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""),comment.char = "#")
  rownames(tRNA_optimal) = tRNA_optimal$codon
  data11[data11$species == species,c("POC1","POC2")] = tRNA_optimal[ data11[data11$species == species,]$codon,c("POC1","POC2")]
}  



data11$color = sapply(data11$codon,function(x)substr(x,3,3))

data11$amino_acid = factor(data11$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa))
data11$anticodon = str_replace_all(data11$anticodon,'T','U')
data11$codon = str_replace_all(data11$codon,'T','U')

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
data11$title = paste(data11$anticodon,"(",data11$codon,")",sep="")
data11$codon = factor(data11$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
data11$title = factor(data11$title,levels= tapply(data11$title, as.integer(data11$codon),unique))

nb_sp = length(unique(data11$species))
data11$nb_species_0 = tapply(data11$tRNA_gene_copy == 0,data11$codon,sum)[data11$codon]
data11$nb_species_0 = round(data11$nb_species_0 / nb_sp*100)
data11$y_axis_0 = tapply(data11$tRNA_gene_copy ,data11$codon,function(x) quantile(x,0.9))[data11$codon]
data11[duplicated(data11$codon) ,]$nb_species_0 = NA
data11[!is.na(data11$nb_species_0),]$nb_species_0 = paste(data11[!is.na(data11$nb_species_0),]$nb_species_0 ,"%")


write.table(data11,"data/data11_supp.tab",quote=F,row.names = F,sep="\t")


