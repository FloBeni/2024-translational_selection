# Generate Data 14
path = "/home/fbenitiere/data/"
path = "/beegfs/data/fbenitiere/"

library(stringi)
stderror <- function(x) sd(x , na.rm = T)/sqrt(length(x[!is.na(x)] ))


code = read.delim(paste(path,"Projet-SplicedVariants/Fichiers-data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")




amino_acid_list = unique(code[code$nb_syn > 1 & code$aa_name != "Ter",]$aa_name_scu)
amino_acid_list = unique(code[ code$aa_name != "Ter",]$aa_name)


data1 = read.delim("data/data1.tab")
data1 = data1[data1$nb_gene > 5000 & data1$pval_aa_fpkm < 0.05,]
list_species = unique(data1$species)


i = 0
data.frame_graph = data.frame()
for( species in list_species ){ print(species)
  i=i+1
  print(i)
  
  codon_usage = read.delim( paste(path,"/Projet-SplicedVariants/Analyses/",species,"/codon_usage_gene_fpkm.tab",sep="") )
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  codon_usage = codon_usage[codon_usage$length_cds == codon_usage$length,]
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage$stop_codon = rowSums(codon_usage[,stop_codon])
  if (!species %in% c("Melospiza_melodia","Oenanthe_oenanthe","Eudyptes_filholi")){codon_usage = codon_usage[codon_usage$stop_codon == 1,]}
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 & !is.na(codon_usage$median_fpkm) ,]
  
  
  if( nrow(codon_usage) == 0){filtering = "not enough CDS"  } else {
    
    ######## EXPECTED CODON OPTIMAL         ########         ########         ########         ########         ########         ########         ######## 
    
    
    for ( type_aa in amino_acid_list){
      codon_used = code[code$aa_name == type_aa,]$codon
      
      aa_per_gene = rowSums(codon_usage[ codon_used ],na.rm = T)
      aa_per_genei = rowSums(codon_usage[ paste(codon_used,'_intronic',sep = "") ],na.rm = T)
      aa_sum = sum(aa_per_gene)
      aa_sumi = sum(aa_per_genei , na.rm=T)
      
      for (codon in codon_used){
        per_gene = codon_usage[,codon]
        per_genei = codon_usage[,paste(codon,'_intronic',sep = "")]
        sum = sum(per_gene) / aa_sum
        sumi = sum(per_genei,na.rm=T) / aa_sumi
        
        data.frame_graph = rbind(data.frame_graph,data.frame(
          species,
          type_aa,
          codon,
          per_gene = mean(per_gene / aa_per_gene, na.rm = T),
          per_genei = mean(per_genei / aa_per_genei , na.rm = T),
          sum,
          sumi
        ))
      }
    }
  }
}

write.table(data.frame_graph,paste("data/data14.tab",sep=""),sep="\t",quote=F,row.names=F)



