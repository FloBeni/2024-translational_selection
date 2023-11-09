# Generate Data 2
library(stringi)
library(ape)

path = "/home/fbenitiere/data/"



code = read.delim(paste(path,"Projet-SplicedVariants/Fichiers-data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )

wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")


data1 = read.delim("data/data1.tab")
data1 = data1[ data1$pval_aa_fpkm < 0.05,]
list_species = unique(data1$species)


all_dt_about_codon = data.frame()
for (species in list_species){print(species)
  if (file.exists(paste(path,"Projet-NeGA/translational_selection/ExpOpti/",species,"_code_table.tab",sep=""))){
    tRNA_optimal = read.delim(paste(path,"Projet-NeGA/translational_selection/ExpOpti/",species,"_code_table.tab",sep=""))
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
}


data2 = data.frame()
for ( codon in unique(code$codon)){
  total = nrow(all_dt_about_codon[all_dt_about_codon$codon == codon,])
  dc = data.frame(table(all_dt_about_codon[all_dt_about_codon$codon == codon,"categorie"]))
  dc$Prop = dc$Freq / total
  dc$total = total
  dc$amino_acid = paste(code[codon,]$aa_name , " (",code[codon,]$nb_syn,")",sep="")
  dc$codon = codon
  dc$WB_type =  wobble_type[substr(codon,3,3)]
  dc$species = "metazoa"
  
  data2 = rbind(data2,dc)
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
    
    data2 = rbind(data2,dc)
  }
}

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
data2$codon = factor(data2$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","T","A","G"),sep=""))) )

data2$title = factor(paste(data2$codon," (",data2$WB_type,")",sep=""),
                     sapply(levels(data2$codon),function(x) paste(x," (",wobble_type[substr(x,3,3)],")",sep="")) )

set_color = c("WCp + abond" = "#33A02C" ,"WCp" = "#B2DF8A","WBp + abond" = "#E31A1C","WBp" = "#FB9A99","not decoded" = "#e2cc1a")
data2$Var1 = factor(data2$Var1,levels =  names(set_color))

write.table(data2,"data/data2.tab",quote=F,row.names = F,sep="\t")





dc = all_dt_about_codon[grepl("abond",all_dt_about_codon$categorie),]
dc$nb_opti = table(paste(dc$aa_name,dc$species,sep="_"))[paste(dc$aa_name,dc$species,sep="_")]
dc$nb_syn = code[dc$codon,]$nb_syn

dc = dc[dc$nb_syn != dc$nb_opti,]

df = data.frame(
  species = names(tapply(dc$codon, dc$species,function(x) sum(substr(x,3,3) %in% c("G","C")) / length(x)  )),
  GC_opti = tapply(dc$codon, dc$species,function(x) sum(substr(x,3,3) %in% c("G","C")) / length(x)  )
)
rownames(df) = df$species




clade_dt = read.delim(paste( "data/clade_dt.tab",sep=""),header=T)
rownames(clade_dt) = clade_dt$species
clade_dt$clade_group = factor(clade_dt$clade_group, levels = c("Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Teleostei","Mammalia","Aves","Other Tetrapodes"))

clade_dt = clade_dt[ clade_dt$clade_group %in% c("Lepido Diptera") & clade_dt$species %in% list_species,]

for (species in clade_dt$species){
  differential_expressed = read.delim(paste(path , "Projet-NeGA/translational_selection/rscu_expressed_genes/intronic_all/" , species , "_DUC" , ".tab" , sep = "" ))
  data_extract = differential_expressed[differential_expressed$COA,]$codon
  GC_DUC =  sum(substr(data_extract , 3 , 3) %in% c("G" , "C")) / length(data_extract)
  df[species , "GC_DUC"] = GC_DUC
}



write.table(df , "data/data13.tab",quote=F,row.names = F,sep="\t")

