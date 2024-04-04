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


wobble_rule = list("G"=c("C","T"),
                   "A"=c("C","T","A"),
                   "T"=c("A","G"),
                   "C"=c("G"))


list_species = read.delim("data/GTDrift_list_species.tab")
rownames(list_species) = list_species$species


for( species in list_species$species ){
  print(species)
  genome_assembly = list_species[species,]$assembly_accession
  taxID = list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  code_table = code
  
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
  
  code_table$nb_tRNA_copies = tRNASE_copies_table[code_table$codon]
  code_table[is.na(code_table$nb_tRNA_copies),]$nb_tRNA_copies = 0
  
  code_table$decoded = F
  code_table$POC1 = F
  code_table$POC2 = F
  
  for (aa in unique(code_table$aa_name)){
    if (length(code_table[code_table$aa_name == aa & code_table$nb_tRNA_copies != 0,]$nb_tRNA_copies) != 0){
      tRNA_present = code_table[code_table$aa_name == aa & code_table$nb_tRNA_copies != 0,]$nb_tRNA_copies
      names(tRNA_present) = code_table[code_table$aa_name == aa & code_table$nb_tRNA_copies != 0,]$anticodon
      
      decoded_codon = unique(unlist(sapply(names(tRNA_present),function(x) paste(chartr( "TUACG" , "AATGC" , stri_reverse(substr(x,2,3))),unlist(wobble_rule[substr(x,1,1)]),sep=""))))
      code_table[code_table$codon %in% decoded_codon,]$decoded = T
      
      if ( all(codon_data[codon_data$aa_name == aa ,]$decoded)){
        if ( !aa %in% c("Ter","Met","Trp")){
          if (length(tRNA_present) > 1 & any(tRNA_present != max(tRNA_present))){
            abundant = names(tRNA_present[ tRNA_present == max(tRNA_present)])
            code_table[code_table$anticodon %in% abundant,]$POC1 = T
            decoded_codon = unique(unlist(sapply(abundant,function(x) paste(chartr( "TUACG" , "AATGC" , stri_reverse(substr(x,2,3))),unlist(wobble_rule[substr(x,1,1)]),sep=""))))
            if (any(code_table$codon %in% decoded_codon & code_table$nb_tRNA_copies == 0)){
              code_table[code_table$codon %in% decoded_codon & code_table$nb_tRNA_copies == 0,]$POC1 = T}
          } else if (length(tRNA_present) == 1){
            abundant = names(tRNA_present[ tRNA_present == max(tRNA_present)])
            code_table[code_table$anticodon %in% abundant,]$POC2 = T
          }
        }
      }
    }}
  print(table(code_table$decoded))
  
  code_table$nb_tRNA_copies = sapply(code_table$codon,function(codon) sum(tRNASE_copies_table[codon],na.rm = T))
  code_table$nb_tRNA_copies_aa = sapply(code_table$aa_name,function(aa) sum(code_table[code_table$aa_name == aa,]$nb_tRNA_copies,na.rm = T))
  code_table$RTF = code_table$nb_tRNA_copies / (code_table$nb_tRNA_copies_aa / code_table$nb_syn)
  
  write.table(code_table,paste(path,"/decoding_table.tab",sep=""),sep="\t",quote=F,row.names=F)
  bash_command <- paste("gzip -f ",path,"/decoding_table.tab",sep="")
  system(bash_command)
}
