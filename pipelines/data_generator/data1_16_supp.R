# Generate Data 1 and 16
library(stringi)


data_conservation = read.delim(paste("data/compilation_prop_gap_pergene_25_50_75.tab.gz",sep=""))
data_conservation_rmfirst1000bp = read.delim(paste("data/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab.gz",sep=""))

stderror <- function(x) sd(x , na.rm = T)/sqrt(length(x[!is.na(x)] ))

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])


wobble_pairing = c("C"="IC","T"="GU","G"="UG","A"="IA")
wobble_associat_wc = list("T"="C", "A"="A", "G"="A", "C"="T")


GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

data1 = data.frame()
data16 = data.frame()
for (species in GTDrift_list_species$species){
  dt_species = data.frame()
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  if (
    file.exists(paste(path,"/trna_from_gff.txt.gz",sep="")) &
    file.size(paste(path,"/trna_from_gff.txt.gz",sep="")) != 38 ){ # if non empty
    tRNA_GFF = T
  } else if (
    file.exists(paste(path,"/trna_from_trnascanse.txt.gz",sep="")) &
    file.size(paste(path,"/trna_from_trnascanse.txt.gz",sep="")) != 36  ){ # if non empty
    tRNA_GFF = F
  } 
  
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.txt.gz",sep="") )
  nb_genes = length(unique(codon_usage$gene_id))
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  
  codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]
  
  if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
    codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}
  
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  prop_cds_expressed = nrow(codon_usage[codon_usage$median_fpkm != 0 ,]) / nrow(codon_usage)
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
  
  observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )
  
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
  
  ## Selectionner les tRNA + abondant pour les duet XXA-XXG
  abond_AG = tRNA_optimal[!tRNA_optimal$aa_name %in% c("Met","Trp","Ile"),]
  abond_AG = abond_AG[substr(abond_AG$codon,3,3) %in% c("A","G"),]
  abond_AG = abond_AG[order(abond_AG$nb_tRNA_copies,decreasing = T),]
  # abond_AG = abond_AG[ !abond_AG$aa_name_scu %in% abond_AG[abond_AG$nb_tRNA_copies == 0,]$aa_name_scu,]
  abond_AG = abond_AG[ abond_AG$aa_name_scu %in% c("Lys","Glu","Leu_4" ,"Gln","Val","Thr","Pro","Leu_2" ,"Ser_4", "Arg_2", "Ala" ),]
  vector = (tapply(abond_AG$nb_tRNA_copies , abond_AG$aa_name_scu,function(x) sum(x == max(x))))
  abond_AG = abond_AG[abond_AG$aa_name_scu %in% names(vector[vector != 2]),]
  abond_AG = abond_AG[!duplicated(paste(abond_AG$aa_name_scu)),]$codon
  print(length(abond_AG))
  
  GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
  ATGCi_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)
  
  GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
  ATGC3_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)
  
  gc3 = (GC3_obs / ATGC3_obs)
  gci = (GCi_obs / ATGCi_obs)
  gci_expressed_10percent = (gci)[codon_usage$median_fpkm >= quantile(codon_usage$median_fpkm , probs = 0.9)]
  model_lm = lm( gc3 ~  gci)
  
  spearman_method_gc3gci = cor.test( gc3, gci,method = "spearman",exact=F)
  
  xaxis = codon_usage$median_fpkm 
  proportion = 2/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  # print(mean(table(intervalle_FPKM)))
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  ## Translational selection signal
  data_optiplus = data.frame()
  for ( amino_acid in unique(tRNA_optimal[ tRNA_optimal$nb_syn >= 2,]$aa_name)){
    codon_used = rownames(tRNA_optimal[tRNA_optimal$aa_name == amino_acid,])
    amino_acid_count = rowSums(codon_usage[codon_used],na.rm = T)
    triplet_intronic = rowSums(codon_usage[paste(codon_used,'_intronic',sep = "")],na.rm = T)
    for ( codon in codon_used ){
      COA_obs =  unlist(codon_usage[codon])
      COA_neg_obs = amino_acid_count
      
      COA_obs_intronic =  unlist(codon_usage[paste(codon,'_intronic',sep = "")])
      COA_neg_obs_intronic = triplet_intronic
      
      data_optiplus = rbind(data_optiplus,data.frame(
        species,
        amino_acid,
        codon,
        aa_name_scu = tRNA_optimal[codon,]$aa_name_scu,
        expressed_overused_background = (round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) -
                                           round( mean( (COA_obs / COA_neg_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T) , 5)) - (
                                             round(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) -
                                               round( mean( (COA_obs_intronic / COA_neg_obs_intronic)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T),5)
                                           )
      ))
    }
  }
  data_optiplus = data_optiplus[order(data_optiplus$expressed_overused_background,decreasing = T),]
  data_optiplus = data_optiplus[order(data_optiplus$amino_acid,decreasing = F),]
  data_optiplus$rank = unlist(tapply(data_optiplus$expressed_overused_background,data_optiplus$amino_acid , function(x) rev(rank(x))))
  ## Faire mieux 
  data_optiplus$nb_tRNA_copies = tRNA_optimal[data_optiplus$codon,]$nb_tRNA_copies
  data16 = rbind(data16,data_optiplus)
  DUC_IC = data_optiplus[substr(data_optiplus$codon,3,3) %in% c("C","T"),]
  DUC_IC = DUC_IC[ DUC_IC$aa_name_scu %in% DUC_IC[DUC_IC$nb_tRNA_copies == 0,]$aa_name_scu,]
  DUC_IC = DUC_IC[!duplicated(DUC_IC$aa_name_scu),]$codon
  table(substr(DUC_IC,3,3))
  
  
  dt_species = rbind(dt_species,data.frame(
    species,
    tRNA_GFF,
    prop_cds_expressed,
    nb_genes,
    nb_genes_filtered = nrow(codon_usage),
    nb_codon_not_decoded,
    
    rho_aa_fpkm = spearman_method_aa$estimate,
    pval_aa_fpkm = spearman_method_aa$p.value,
    
    rho_gc3_gci = spearman_method_gc3gci$estimate,
    pval_gc3_gci = spearman_method_gc3gci$p.value,
    
    g_abond_ag = sum(substr(abond_AG,3,3) %in% c("G")) / length(abond_AG),
    c_duc_ic = sum(substr(DUC_IC,3,3) %in% c("C")) / length(DUC_IC),
    
    gc3 = mean(gc3,na.rm=T),
    var_gc3 = var(gc3,na.rm=T),
    gci = mean(gci,na.rm=T),
    std_gci = stderror(gci),
    var_gci = var(gci,na.rm=T),
    var_gci_exp = var(gci_expressed_10percent,na.rm=T)
    
  ) )
  
  ### Translational selection intensity
  
  for ( type_aa in c( "POC1","POC2","POCs")){
    if (type_aa == "POC1" ){
      subset_selected = tRNA_optimal[tRNA_optimal$POC1,]
      list_POC = subset_selected$codon
      list_aa = unique(subset_selected$aa_name)
      list_codon = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa,]$codon
    } else if  (type_aa == "POC2" ){
      subset_selected = tRNA_optimal[tRNA_optimal$POC2,]
      list_POC = subset_selected$codon
      list_aa = unique(subset_selected$aa_name)
      list_codon = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa,]$codon
    } else if  (type_aa == "POCs" ){
      subset_selected = tRNA_optimal[tRNA_optimal$POC2 | tRNA_optimal$POC1,]
      list_POC = subset_selected$codon
      list_aa = unique(subset_selected$aa_name)
      list_codon = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa,]$codon
    }
    
    
    if ( length(list_POC) != 0 ){
      
      ##### Over-used of POC in expressed genes
      POC_obs = rowSums(codon_usage[ list_POC ],na.rm = T)
      POC_codons_obs = rowSums(codon_usage[ list_codon ],na.rm = T)
      POC_obs_intronic = rowSums(codon_usage[paste(list_POC,'_intronic',sep = "")],na.rm = T)
      POC_codons_obs_intronic = rowSums(codon_usage[paste(list_codon,'_intronic',sep = "")],na.rm = T)
      
      
      ## Over-used of POC in constraint sites
      
      if (GTDrift_list_species[species,]$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei")){
        data_conservation_sub = data_conservation_rmfirst1000bp[data_conservation_rmfirst1000bp$species == species & data_conservation_rmfirst1000bp$protein %in% codon_usage$protein_id,]
      } else {
        data_conservation_sub = data_conservation[data_conservation$species == species ,]
      }
      
      
      table_constrain = data.frame(busco_id = data_conservation_sub$busco_id)
      for (constrain in c("_highconst","_modconst","_sligconst","_unconst")){
        table_constrain[,paste("POC",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_POC,constrain,sep = "")],na.rm = T)
        table_constrain[,paste("POC_codon",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_codon,constrain,sep = "")],na.rm = T)
      }
      
      poc_exp_genes = (POC_obs / POC_codons_obs)[intervalle_FPKM == names(FPKM_bins)[length(FPKM_bins)]]
      poc_noexp_genes = (POC_obs / POC_codons_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )]
      
      Fpoc_expressed = round(mean(poc_exp_genes,na.rm = T),5)
      Fpoc_noexpressed = round(mean(poc_noexp_genes,na.rm = T),5)
      
      dt_inter_err = data.frame()
      for (i in 1:1000){
        boot_poc_exp_genes =  sample(poc_exp_genes, replace = T)
        boot_poc_noexp_genes =  sample(poc_noexp_genes, replace = T)
        
        dt_inter_err = rbind(dt_inter_err,data.frame(
          S = log( mean( sample(boot_poc_exp_genes, replace = T),na.rm=T)/(1 - mean( sample(boot_poc_exp_genes, replace = T),na.rm=T))) - log(mean( sample(boot_poc_noexp_genes, replace = T),na.rm=T)/(1 - mean( sample(boot_poc_noexp_genes, replace = T),na.rm=T)))
        ))
      }
      
      
      dt_translational_selection = data.frame(
        nb_aa = length(list_aa),
        nb_genes_per_bins = mean(table(intervalle_FPKM)),
        nb_busco = nrow(table_constrain),
        S = log(Fpoc_expressed/(1-Fpoc_expressed)) - log(Fpoc_noexpressed/(1-Fpoc_noexpressed)),
        S_int_025 = quantile(dt_inter_err$S,c(0.025)),
        S_int_975 = quantile(dt_inter_err$S,c(0.975)),
        expressed_overused = 100 * (Fpoc_expressed - Fpoc_noexpressed) ,
        expressed_overused_background = 100 * ((Fpoc_expressed - Fpoc_noexpressed) - (
          round(tapply( POC_obs_intronic / POC_codons_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) -
            round( mean( (POC_obs_intronic / POC_codons_obs_intronic)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T),5)
        )),
        # constraint_overused = 100*(mean(table_constrain$POC_highconst/table_constrain$POC_codon_highconst ,na.rm = T) - mean(table_constrain$POC_unconst/table_constrain$POC_codon_unconst,na.rm = T))
        constraint_overused = 100 * mean(table_constrain$POC_highconst/table_constrain$POC_codon_highconst - table_constrain$POC_unconst/table_constrain$POC_codon_unconst,na.rm = T)
      )
    } else {
      dt_translational_selection = data.frame(
        nb_aa = 0,
        nb_genes_per_bins = NA,
        nb_busco = NA,
        S = NA,
        S_int_025 = NA,
        S_int_975 = NA,
        expressed_overused = NA,
        expressed_overused_background = NA,
        constraint_overused = NA
      )
    }
    colnames(dt_translational_selection) = paste(colnames(dt_translational_selection),type_aa,sep="_")
    
    dt_species = cbind(dt_species,dt_translational_selection)
  }
  
  data1 = rbind(data1,dt_species)
}
write.table(data1,"data/data1_supp.tab",quote=F,row.names = F,sep="\t")

data16 = data16[data16$species %in% GTDrift_list_species[GTDrift_list_species$clade_group %in% c("Diptera","Lepidoptera"),]$species,]
write.table(data16,"data/data16.tab",quote=F,row.names = F,sep="\t")

