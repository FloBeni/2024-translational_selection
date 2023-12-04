# Download from source
options(stringsAsFactors = F, scipen = 999)


pathData="/home/fbenitiere/data/Projet-SplicedVariants/"
pathData="/beegfs/data/fbenitiere/Projet-SplicedVariants/"

bash_command <- paste("mkdir -p ","data/per_species/",sep="")
system(bash_command)


list_species = read.delim("data/GTDrift_list_species.tab")
rownames(list_species) = list_species$species

for (species in list_species$species){print(species)
  genome_assembly = list_species[species,]$assembly_accession
  taxID = list_species[species,]$NCBI.taxid
  
  bash_command <- paste("mkdir -p data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  system(bash_command)
  
  
  # Get tRNA table and compress.
  bash_command <- paste("cp ",pathData,"Annotations/",species,"/formatted_data/tRNAscan.tab ","data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/tRNA_from_GFF.tab",sep="")
  system(bash_command)
  bash_command <- paste("gzip data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/tRNA_from_GFF.tab",sep="")
  system(bash_command)
  
  if (file.size(paste(pathData,"Annotations/",species,"/formatted_data/tRNAscan_SE.tab",sep="")) != 0){
    tRNASE_copies = read.delim(paste(pathData,"Annotations/",species,"/formatted_data/tRNAscan_SE.tab",sep=""), header = F)
    colnames(tRNASE_copies) = unlist(tRNASE_copies[2,])
    tRNASE_copies = tRNASE_copies[ !is.na(as.numeric(tRNASE_copies$Score)),]
    tRNASE_copies = tRNASE_copies[,-c(2)]
    write.table(tRNASE_copies,paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/tRNAscan_SE.tab",sep=""),sep="\t",quote=F,row.names=F)
  } else {
    bash_command <- paste("cp ",pathData,"Annotations/",species,"/formatted_data/tRNAscan_SE.tab ","data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
    system(bash_command)
  }
  bash_command <- paste("gzip data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/tRNAscan_SE.tab",sep="")
  system(bash_command)
  
  
  bash_command <- paste("cp ",pathData,"Analyses/",species,"/codon_usage_gene_fpkm.tab ","data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  system(bash_command)
  bash_command <- paste("gzip data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/codon_usage_gene_fpkm.tab",sep="")
  system(bash_command)
}

