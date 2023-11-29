# Generate Data 2
library(stringi)
library(ape)


code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )

wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


data1 = read.delim("data/data1.tab")
data1 = data1[ data1$pval_aa_fpkm < 0.05 ,]


all_dt_about_codon = data.frame()
for (species in unique(data1$species)){
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  rownames(tRNA_optimal) = tRNA_optimal$codon
  
  tRNA_optimal[tRNA_optimal$nb_tRNA_copies != 0,"categorie"] = "WCp"
  tRNA_optimal[tRNA_optimal$nb_tRNA_copies == 0,"categorie"] = "WBp"
  tRNA_optimal[tRNA_optimal$WC_abond,"categorie"] = "WCp + abond"
  tRNA_optimal[tRNA_optimal$Wobble_abond,"categorie"] = "WBp + abond"
  tRNA_optimal[!tRNA_optimal$decoded,"categorie"] = "not decoded"
  
  tRNA_optimal = tRNA_optimal[,c("codon","aa_name","anticodon","categorie")]
  tRNA_optimal$species = species
  
  all_dt_about_codon = rbind(all_dt_about_codon,tRNA_optimal)
}


data4 = data.frame()
for ( codon in unique(code$codon)){
  total = nrow(all_dt_about_codon[all_dt_about_codon$codon == codon,])
  dc = data.frame(table(all_dt_about_codon[all_dt_about_codon$codon == codon,"categorie"]))
  dc$Prop = dc$Freq / total
  dc$total = total
  dc$amino_acid = paste(code[codon,]$aa_name , " (",code[codon,]$nb_syn,")",sep="")
  dc$codon = codon
  dc$WB_type =  wobble_type[substr(codon,3,3)]
  dc$species = "metazoa"
  
  data4 = rbind(data4,dc)
}


for (species in c("Homo_sapiens","Caenorhabditis_elegans","Drosophila_melanogaster")){
  all_dt_about_codon_sub = all_dt_about_codon[all_dt_about_codon$species == species  ,]
  
  for ( codon in unique(code$codon)){
    total = nrow(all_dt_about_codon_sub[all_dt_about_codon_sub$codon == codon,])
    dc = data.frame(table(all_dt_about_codon_sub[all_dt_about_codon_sub$codon == codon,"categorie"]))
    dc$Prop = dc$Freq / total
    dc$total = total
    dc$amino_acid = paste(code[codon,]$aa_name , " (",code[codon,]$nb_syn,")",sep="")
    dc$codon = codon
    dc$WB_type =  wobble_type[substr(codon,3,3)]
    dc$species = species
    
    data4 = rbind(data4,dc)
  }
}

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
data4$codon = factor(data4$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","T","A","G"),sep=""))) )

data4$title = factor(paste(data4$codon," (",data4$WB_type,")",sep=""),
                     sapply(levels(data4$codon),function(x) paste(x," (",wobble_type[substr(x,3,3)],")",sep="")) )

set_color = c("WCp + abond" = "#33A02C" ,"WCp" = "#B2DF8A","WBp + abond" = "#E31A1C","WBp" = "#FB9A99","not decoded" = "#e2cc1a")
data4$Var1 = factor(data4$Var1,levels =  names(set_color))

write.table(data4,"data/data4.tab",quote=F,row.names = F,sep="\t")






