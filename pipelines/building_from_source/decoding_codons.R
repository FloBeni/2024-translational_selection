options( stringsAsFactors = F, scipen = 999 )
library(stringi)

path = "/home/fbenitiere/data/"
# path = "/beegfs/data/fbenitiere/"

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")


wobble_rule = list("T"=c("A","G"),
                   "A"=c("T","A"),
                   "G"=c("C","T"),
                   "C"=c("G","A"))


list_species = read.delim("data/GTDrift_list_species.tab")
rownames(list_species) = list_species$species


for( species in list_species$species ){
  print(species)
  genome_assembly = list_species[species,]$assembly_accession
  taxID = list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  if (
    file.exists(paste(path,"/tRNA_from_GFF.tab.gz",sep="")) &
    file.size(paste(path,"/tRNA_from_GFF.tab.gz",sep="")) != 38 ){
    tRNASE_gff = read.delim( paste(path,"/tRNA_from_GFF.tab.gz",sep="") )
    tRNASE_gff$codon = sapply(tRNASE_gff$anticodon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
    tRNASE_gff_table = table(tRNASE_gff$codon)
    tRNASE_copies_table = tRNASE_gff_table
    tRNA_GFF = T
  } else if (
    file.exists(paste(path,"/tRNAscan_SE.tab.gz",sep="")) &
    file.size(paste(path,"/tRNAscan_SE.tab.gz",sep="")) != 36  ){
    tRNASE_copies = read.delim(paste(path,"/tRNAscan_SE.tab.gz",sep=""), header = T)
    tRNASE_copies = tRNASE_copies[as.numeric(tRNASE_copies$Score) > 55,]
    tRNASE_copies = tRNASE_copies[tRNASE_copies$Note != "pseudo" | is.na(tRNASE_copies$Note),]
    tRNASE_copies$anticodon = sapply(tRNASE_copies$Codon,function(x) chartr("TUACG","TTACG",x)  )
    tRNASE_copies$codon = sapply(tRNASE_copies$Codon,function(x) chartr("TUACG","AATGC",stri_reverse(x)) )
    tRNASE_copies_table = table(tRNASE_copies$codon)
    tRNA_GFF = F
  } 
  
  
  codon_data = data.frame()
  for (codon in unique(code$codon)){
    codon_data = rbind(codon_data,
                       data.frame(
                         codon = codon,
                         anticodon = chartr("TUACG","AATGC", stri_reverse(codon)),
                         tRNA_WC_copies= sum(tRNASE_copies_table[codon],na.rm = T)
                       ))
  }
  codon_data$amino_acid = code[codon_data$codon,]$aa_name
  
  codon_data$letter_aa = code[codon_data$codon,]$aa
  codon_data = codon_data[order(codon_data$tRNA_WC_copies,decreasing = T),]
  
  # Identify abundant tRNA
  aa_data = data.frame()
  for ( amino_acid in unique(code$aa_name)){
    codon_selected = codon_data[codon_data$amino_acid == amino_acid, ]
    tRNA_available_nb = nrow(codon_selected)
    antic_abondant = codon_selected[codon_selected$tRNA_WC_copies == max(codon_selected$tRNA_WC_copies) & 
                                      codon_selected$tRNA_WC_copies != 0  , ]
    
    if ( nrow(antic_abondant) == 0 | nrow(antic_abondant) == tRNA_available_nb ){ 
      # print("pas de tRNA plus abondant que les autres")
      # print(amino_acid)
    } else {
      aa_data = rbind(aa_data,
                      data.frame(
                        amino_acid ,
                        letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                        tRNA_copies= antic_abondant$tRNA_WC_copies,
                        codon =  antic_abondant$codon,
                        anticodon =  antic_abondant$anticodon
                      ))
    }
  }
  list_OExp_WC = aa_data$codon
  code_table = code
  code_table$WC_abond = F
  code_table[list_OExp_WC,]$WC_abond = T
  
  
  # Identify codon wobble decoded by an abundant tRNA
  list_OExp_Wobble = c()
  for ( amino_acid in aa_data$amino_acid){
    codon_used = rownames(code[code$aa_name == amino_acid,])
    antic_abondant = aa_data[aa_data$amino_acid  == amino_acid,]$anticodon
    for ( codon in codon_used ){
      anticodon = paste(unlist(wobble_rule[substr(codon,3,3)]),
                        chartr( "TUACG" , "AATGC" , stri_reverse(substr(codon,1,2))),sep = "")
      if ( codon_data[codon_data$codon == codon ,]$tRNA_WC_copies == 0 &
           any( anticodon %in% antic_abondant ) ){
        list_OExp_Wobble = append(list_OExp_Wobble,codon)
      }
    }
  }
  
  code_table$Wobble_abond = F
  if (length(list_OExp_Wobble) != 0){
    code_table[list_OExp_Wobble,]$Wobble_abond = T
  }
  
  # Identify codon decoded by any tRNA
  list_Exp = c()
  for ( amino_acid in unique(code$aa_name)){
    codon_used = rownames(code[code$aa_name == amino_acid,])
    antic_abondant = codon_data[codon_data$amino_acid  == amino_acid & codon_data$tRNA_WC_copies != 0,]$anticodon
    for ( codon in codon_used ){
      anticodon = paste(unlist(wobble_rule[substr(codon,3,3)]),
                        chartr( "TUACG" , "AATGC" , stri_reverse(substr(codon,1,2))),sep = "")
      if ( codon_data[codon_data$codon == codon ,]$tRNA_WC_copies == 0 &
           any( anticodon %in% antic_abondant ) ){
        list_Exp = append(list_Exp,codon)
      } else if ( codon_data[codon_data$codon == codon ,]$tRNA_WC_copies != 0 ){
        list_Exp = append(list_Exp,codon)
      }
    }
  }
  
  code_table$decoded = F
  if (length(list_Exp) != 0){
    code_table[list_Exp,]$decoded = T
  }
  
  
  code_table$nb_tRNA_copies = sapply(code_table$codon,function(codon) sum(tRNASE_copies_table[codon],na.rm = T))
  code_table$nb_tRNA_copies_aa = sapply(code_table$aa_name,function(aa) sum(codon_data[codon_data$amino_acid == aa,]$tRNA_WC_copies,na.rm = T))
  code_table$RTF = code_table$nb_tRNA_copies / (code_table$nb_tRNA_copies_aa / code_table$nb_syn)
  rownames(code_table) = code_table$anticodon
  
  ## POC1 : Codons decoded in Watson-Crick or Wobble of the most abundant tRNA(s).
  nb_decoded_most_abundant = table(code_table[(code_table$WC_abond | code_table$Wobble_abond), ]$aa_name)
  aa_POC1 = names(nb_decoded_most_abundant[nb_decoded_most_abundant != table(code_table$aa_name)[names(nb_decoded_most_abundant)]])
  code_table$POC1 = (code_table$WC_abond | code_table$Wobble_abond) & code_table$aa_name %in% aa_POC1
  
  ## POC2 : Codons decoded in Watson-Crick of the unique tRNA available.
  nb_trna_0 = table(code_table[code_table$nb_tRNA_copies != 0, ]$aa_name)
  aa_POC2 = names(nb_trna_0)[nb_trna_0 == 1]
  code_table$POC2 = code_table$WC_abond & code_table$aa_name %in% aa_POC2
  
  print(table(code_table$POC1 & code_table$POC2))
  
  write.table(code_table,paste(path,"/decoding_table.tab",sep=""),sep="\t",quote=F,row.names=F)
  bash_command <- paste("gzip -f ",path,"/decoding_table.tab",sep="")
  system(bash_command)
}



# unique(data_by_species[!data_by_species$decoded & data_by_species$aa_name != "Ter",]$species)
# 
# data_by_species_decoded = data_by_species[!data_by_species$species %in% data_by_species[!data_by_species$decoded & data_by_species$aa_name != "Ter",]$species,]
# 
# mean(table(data_by_species_decoded[data_by_species_decoded$POC1,]$species))
# mean(table(data_by_species_decoded[data_by_species_decoded$POC2,]$species))
# table(data_by_species_decoded[data_by_species_decoded$POC2,]$codon)
# 
# table(data_by_species_decoded$POC1 & data_by_species_decoded$POC2 & !data_by_species_decoded$aa_name=="Ter")

