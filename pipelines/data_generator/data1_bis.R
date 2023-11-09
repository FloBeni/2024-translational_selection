# Generate Data 1
library(stringi)

stderror <- function(x) sd(x , na.rm = T)/sqrt(length(x[!is.na(x)] ))

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


data1 = data.frame()
for (species in GTDrift_list_species$species){
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxID",taxID,"/",genome_assembly,sep="")
  
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )
  nb_genes = nrow(codon_usage)
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  
  codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]
  
  if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
    codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}
  
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
  
  observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )
  
  GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
  ATGC3_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)
  
  GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
  ATGCi_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)
  
  gc3 = (GC3_obs / ATGC3_obs)
  gci = (GCi_obs / ATGCi_obs)
  model_lm = lm( gc3 ~  gci)
  
  spearman_method_gc3gci = cor.test( gc3, gci,method = "spearman",exact=F)
  
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  rownames(tRNA_optimal) = tRNA_optimal$codon
  tRNA_optimal = tRNA_optimal[tRNA_optimal$aa_name != "Ter",]
  nb_codon_not_decoded = sum(!tRNA_optimal$decoded)
  
  
  aa_data = data.frame()
  for (amino_acid in unique(code$aa_name)){
    codon_used = rownames(code[code$aa_name == amino_acid,])
    aa_data = rbind(aa_data,
                    data.frame(
                      amino_acid ,
                      letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                      tRNASE_copies= sum(tRNA_optimal[codon_used,]$nb_tRNA_copies,na.rm = T),
                      obs_codon = sum(observation[codon_used])
                    ))
  }
  aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]
  
  spearman_method_aa = cor.test( aa_data$tRNASE_copies, aa_data$obs_codon,method="spearman",exact=F)
  
  
  data1 = rbind(data1,data.frame(
    species,
    nb_genes,
    nb_genes_filtered = nrow(codon_usage),
    nb_codon_not_decoded,
    
    rho_aa_fpkm = spearman_method_aa$estimate,
    pval_aa_fpkm = spearman_method_aa$p.value,
    
    
    rho_gc3_gci = spearman_method_gc3gci$estimate,
    pval_gc3_gci = spearman_method_gc3gci$p.value,
    
    gc3 = mean(GC3_obs / ATGC3_obs,na.rm=T),
    std_gc3 = stderror(GC3_obs / ATGC3_obs),
    var_gc3 = var(GC3_obs / ATGC3_obs,na.rm=T),
    
    gci = mean(GCi_obs / ATGCi_obs,na.rm=T),
    std_gci = stderror(GCi_obs / ATGCi_obs),
    var_gci = var(GCi_obs / ATGCi_obs,na.rm=T)
    
    
  ))
  
}

write.table(data1,"data/data1_bis.tab",quote=F,row.names = F,sep="\t")