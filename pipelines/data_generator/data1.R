# Generate Data 1
library(stringi)

path = "/home/fbenitiere/data/"

code = read.delim(paste(path,"Projet-SplicedVariants/Fichiers-data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")


clade_dt = read.delim(paste( "data/clade_dt.tab",sep=""),header=T)
list_species = clade_dt$species

data1 = data.frame()
for (species in list_species){
  # print(species)
  codon_usage = read.delim( paste(path,"/Projet-SplicedVariants/Analyses/",species,"/codon_usage_gene_fpkm.tab",sep="") )
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  
  codon_usage = codon_usage[codon_usage$intern_stop_codon==0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %%3 == 0,]
  
  if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){ 
    codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else {    print(species)}
  
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  
  if (
    file.exists(paste(path,"/Projet-SplicedVariants/Annotations/",species,"/formatted_data/tRNAscan.tab",sep="")) &
    file.size(paste(path,"/Projet-SplicedVariants/Annotations/",species,"/formatted_data/tRNAscan.tab",sep="")) != 0 ){
    tRNASE_gff = read.delim( paste(path,"/Projet-SplicedVariants/Annotations/",species,"/formatted_data/tRNAscan.tab",sep="") )
    tRNASE_gff$codon = sapply(tRNASE_gff$anticodon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
    tRNASE_gff_table = table(tRNASE_gff$codon)
    tRNASE_copies_table = tRNASE_gff_table
    tRNA_GFF = "from GFF"
  } else if (
    file.exists(paste(path,"/Projet-SplicedVariants/Annotations/",species,"/formatted_data/tRNAscan_SE.tab",sep="")) &
    file.size(paste(path,"/Projet-SplicedVariants/Annotations/",species,"/formatted_data/tRNAscan_SE.tab",sep="")) != 0  ){
    tRNASE_copies = read.delim(paste(path,"/Projet-SplicedVariants/Annotations/",species,"/formatted_data/tRNAscan_SE.tab",sep=""), header = F)
    colnames(tRNASE_copies) = unlist(tRNASE_copies[2,])
    tRNASE_copies = tRNASE_copies[4:nrow(tRNASE_copies),]
    tRNASE_copies = tRNASE_copies[as.vector(tRNASE_copies$Note) != "pseudo",]
    tRNASE_copies = tRNASE_copies[tRNASE_copies$Note != "pseudo",]
    tRNASE_copies = tRNASE_copies[as.numeric(tRNASE_copies$Score) > 55 & !is.na(as.numeric(tRNASE_copies$Score)),]
    tRNASE_copies$anticodon = sapply(tRNASE_copies$Codon, function(x) chartr("TUACG","TTACG",x)  )
    tRNASE_copies$codon = sapply(tRNASE_copies$Codon, function(x) chartr("TUACG","AATGC",stri_reverse(x)) )
    tRNASE_copies_table = table(tRNASE_copies$codon)
    tRNA_GFF = "not from GFF"
  } else { tRNASE_copies_table = 0
  tRNA_GFF = "" }
  
  observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )
  
  aa_data = data.frame()
  for (amino_acid in unique(code$aa_name)){
    codon_used = rownames(code[code$aa_name == amino_acid,])
    aa_data = rbind(aa_data,
                    data.frame(
                      species,
                      tRNA_GFF,
                      amino_acid ,
                      letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                      tRNASE_copies= sum(tRNASE_copies_table[codon_used],na.rm = T),
                      obs_codon = sum(observation[codon_used]),
                      nb_gene = nrow(codon_usage)
                    ))
  }
  
  aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]
  spearman_method_aa = cor.test( aa_data$tRNASE_copies, aa_data$obs_codon,method="spearman",exact=F)
  
  aa_data$rho_aa_fpkm = spearman_method_aa$estimate
  aa_data$pval_aa_fpkm = spearman_method_aa$p.value
  data1 = rbind( data1 , aa_data )
}


write.table(data1,"data/data1.tab",quote=F,row.names = F,sep="\t")


