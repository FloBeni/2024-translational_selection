# Generate Data 1
library(stringi)

stderror <- function(x) sd(x , na.rm = T)/sqrt(length(x[!is.na(x)] ))

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])


wobble_pairing = c("C"="IC","T"="GU","G"="UG","A"="IA")
wobble_associat_wc = list("T"="C",
                          "A"="A",
                          "G"="A",
                          "C"="T")


GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


data1 = data.frame()
for (species in GTDrift_list_species$species){
  dt_species = data.frame()
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
  
  
  dt_species = rbind(dt_species,data.frame(
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
    
  ) )
  
  ### Translational selection intensity
  
  tRNA_optimal$decoded_codon = sapply(tRNA_optimal$codon,function(x) paste(substr(x,1,2),wobble_associat_wc[[substr(x,3,3)]],collapse="",sep=""))
  tRNA_optimal$wb_type = wobble_pairing[substr(tRNA_optimal$codon,3,3)]
  
  xaxis = codon_usage$median_fpkm 
  proportion = 5/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  for ( type_aa in c( "WB_WC_notambiguous","WC_duet_ambiguous")){
    if (type_aa == "WB_WC_notambiguous" ){
      dt_selected = tRNA_optimal[ tRNA_optimal$nb_syn >= 2,]
      optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
      
      list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
      nb_aa = length(list_aa)
      
      list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon
      list_COA_neg = code[code$aa_name %in% list_aa,]$codon
      nb_codon = length(list_COA)
      
    }  else if  (type_aa == "WC_duet_ambiguous" ){
      dt_selected = tRNA_optimal[  tRNA_optimal$nb_syn  == 2 ,]
      optimal_count = table( dt_selected[dt_selected$Wobble_abond,]$aa_name )
      
      list_aa = names(optimal_count)[optimal_count != table(code$aa_name)[names(optimal_count)]]
      nb_aa = length(list_aa)
      
      list_COA = dt_selected[dt_selected$WC_abond & dt_selected$aa_name %in% list_aa,]$codon
      list_COA_neg = code[code$aa_name %in% list_aa ,]$codon
      
    }
    
    COA_obs = rowSums(codon_usage[ list_COA ],na.rm = T)
    COA_neg_obs = rowSums(codon_usage[ list_COA_neg ],na.rm = T)
    
    
    if (length(list_COA) != 0){
      COA_obs_intronic = rowSums(codon_usage[paste(list_COA,'_intronic',sep = "")],na.rm = T)
      COA_neg_obs_intronic = rowSums(codon_usage[paste(list_COA_neg,'_intronic',sep = "")],na.rm = T)
    } else {
      COA_obs_intronic = rowSums(codon_usage[list_COA],na.rm = T)
      COA_neg_obs_intronic = rowSums(codon_usage[list_COA_neg],na.rm = T)
    } 
    
    dt_translational_selection = data.frame(
      average_RTF_abundant_tRNA = mean( tRNA_optimal[list_COA,]$RTF[tRNA_optimal[list_COA,]$RTF != 0] ) ,
      expressed_overused = round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) - 
        round( mean( (COA_obs / COA_neg_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T) , 5) ,
      expressed_overused_background = (round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) - 
                                         round( mean( (COA_obs / COA_neg_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T) , 5)) - (
                                                          round(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) -
                                                            round( mean( (COA_obs_intronic / COA_neg_obs_intronic)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T),5)
                                                        ))
    colnames(dt_translational_selection) = paste(colnames(dt_translational_selection),type_aa,sep="_")
    
    dt_species = cbind(dt_species,dt_translational_selection)
  }
  
  data1 = rbind(data1,dt_species)
}
write.table(data1,"data/data1_bis.tab",quote=F,row.names = F,sep="\t")
