# Download from source data on the computing facilities of the CC LBBE/PRABI and the Core Cluster of the Institut Français de Bioinformatique (IFB)
options(stringsAsFactors = F, scipen = 999)


pathData="/home/fbenitiere/gtdrift/data/genome_assembly/"
# pathData="/beegfs/banque/gtdrift/data/genome_assembly/"

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
  bash_command <- paste("cp ",pathData,genome_assembly,"/analyses/codon_usage/trna_pool/trna_from_gff.txt ","data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  system(bash_command)
  bash_command <- paste("gzip data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/trna_from_gff.txt",sep="")
  system(bash_command)
  
  if (file.size(paste(pathData,genome_assembly,"/analyses/codon_usage/trna_pool/trna_from_trnascanse.txt",sep="")) != 0){
    tRNASE_copies = read.delim(paste(pathData,genome_assembly,"/analyses/codon_usage/trna_pool/trna_from_trnascanse.txt",sep=""), header = F)
    colnames(tRNASE_copies) = unlist(tRNASE_copies[2,])
    tRNASE_copies = tRNASE_copies[ !is.na(as.numeric(tRNASE_copies$Score)),]
    tRNASE_copies = tRNASE_copies[,-c(2)]
    write.table(tRNASE_copies,paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/trna_from_trnascanse.txt",sep=""),sep="\t",quote=F,row.names=F)
  } else {
    bash_command <- paste("cp ",pathData,genome_assembly,"/analyses/codon_usage/trna_pool/trna_from_trnascanse.txt ","data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
    system(bash_command)
  }
  bash_command <- paste("gzip data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/trna_from_trnascanse.txt",sep="")
  system(bash_command)
  
  
  bash_command <- paste("cp ",pathData,genome_assembly,"/analyses/codon_usage/codon_usage_gene_fpkm.txt ","data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  system(bash_command)
  bash_command <- paste("gzip data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,"/codon_usage_gene_fpkm.txt",sep="")
  system(bash_command)
}


# Collect data of site constraint 
bash_command <- paste("cp /home/fbenitiere/data/Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75.tab data/compilation_prop_gap_pergene_25_50_75.tab; gzip data/compilation_prop_gap_pergene_25_50_75.tab",sep="")
system(bash_command)

bash_command <- paste("cp /home/fbenitiere/data/Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab data/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab; gzip data/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab",sep="")
system(bash_command)
