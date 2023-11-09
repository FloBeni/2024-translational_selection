# Generate Data 3
library(stringi)
library(ape)


code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")
stop_codon = rownames(code[code$aa_name == "Ter",])

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

### HOMO SAPIENS
species = "Homo_sapiens"
print(species)
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid

path = paste("data/per_species/",species,"_NCBI.taxID",taxID,"/",genome_assembly,sep="")

codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )

codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)

codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]

if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
  codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}

codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]



GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
ATGC3_neg_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)

GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
ATGCi_neg_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)

GC3 = GC3_obs / ATGC3_neg_obs
GCi = GCi_obs / ATGCi_neg_obs

data3 = data.frame()
data3 = rbind(data3,data.frame(
  species,
  GCi,
  GC3,
  categorie = "GC3"
))


#### Caenorhabditis_elegans
species = "Caenorhabditis_elegans"
print(species)
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid

path = paste("data/per_species/",species,"_NCBI.taxID",taxID,"/",genome_assembly,sep="")

codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )

codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)

codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]

if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
  codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}

codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]



GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
ATGC3_neg_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)

GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
ATGCi_neg_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)

GC3 = GC3_obs / ATGC3_neg_obs
GCi = GCi_obs / ATGCi_neg_obs

data3 = rbind(data3,data.frame(
  species,
  GCi,
  GC3,
  categorie = "GC3"
))


#### Drosophila_melanogaster
species = "Drosophila_melanogaster"
print(species)
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid

path = paste("data/per_species/",species,"_NCBI.taxID",taxID,"/",genome_assembly,sep="")

codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )

codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)

codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]

if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
  codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}

codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]



GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
ATGC3_neg_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)

GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
ATGCi_neg_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)

GC3 = GC3_obs / ATGC3_neg_obs
GCi = GCi_obs / ATGCi_neg_obs

data3 = rbind(data3,data.frame(
  species,
  GCi,
  GC3,
  categorie = "GC3"
))


write.table(data3,"data/data3.tab",quote=F,row.names = F,sep="\t")
