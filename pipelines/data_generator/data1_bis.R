# Generate Data 1
library(stringi)

path_data = "/home/fbenitiere/data/"

data_conservation = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75.tab",sep=""))
data_conservation_rmfirst1000bp = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab",sep=""))
# data_conservation = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_cons_pergene_25_50_75.tab",sep=""))
# data_conservation_rmfirst1000bp = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_cons_pergene_25_50_75_rmfirst1000bp.tab",sep=""))

stderror <- function(x) sd(x , na.rm = T)/sqrt(length(x[!is.na(x)] ))

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])
code$nb_syn = table(code$aa_name)[code$aa_name]


wobble_pairing = c("C"="IC","T"="GU","G"="UG","A"="IA")
wobble_associat_wc = list("T"="C",
                          "A"="A",
                          "G"="A",
                          "C"="T")


GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species
# GTDrift_list_species = GTDrift_list_species[GTDrift_list_species$clade_group == "Lepido Diptera",]


data1 = data.frame()
for (species in GTDrift_list_species$species){
  dt_species = data.frame()
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )
  nb_genes = length(unique(codon_usage$gene_id))
  
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
  abond_AG = tRNA_optimal[!tRNA_optimal$aa_name %in% c("Phe","Asn","Asp","His","Cys","Tyr","Met","Trp","Ile"),]
  abond_AG = abond_AG[!substr(abond_AG$codon,3,3) %in% c("C","T"),]
  abond_AG = abond_AG[order(abond_AG$nb_tRNA_copies,decreasing = T),]
  abond_AG = abond_AG[ !abond_AG$aa_name_scu %in% abond_AG[abond_AG$nb_tRNA_copies == 0,]$aa_name_scu,]
  vector = (tapply(abond_AG$nb_tRNA_copies , abond_AG$aa_name_scu,function(x) sum(x == max(x))))
  abond_AG = abond_AG[abond_AG$aa_name_scu %in% names(vector[vector != 2]),]
  abond_AG = abond_AG[!duplicated(paste(abond_AG$aa_name_scu)),]$codon
  
  abond_duet = tRNA_optimal[tRNA_optimal$nb_syn_scu == 2 & tRNA_optimal$aa_name %in% c("Lys" , "Glu" , "Gln" , "Arg" , "Leu"),]
  abond_duet = abond_duet[ !abond_duet$aa_name %in% abond_duet[abond_duet$nb_tRNA_copies == 0,]$aa_name,]
  abond_duet = abond_duet[order(abond_duet$nb_tRNA_copies,decreasing = T),]
  vector = (tapply(abond_duet$nb_tRNA_copies , abond_duet$aa_name,function(x) sum(x == max(x))))
  abond_duet = abond_duet[abond_duet$aa_name %in% names(vector[vector != 2]),]
  abond_duet = abond_duet[!duplicated(abond_duet$aa_name),]$codon
  
  amino_acid_gc2 = tRNA_optimal[tRNA_optimal$nb_syn == 2 & !tRNA_optimal$aa_name %in% tRNA_optimal[tRNA_optimal$Wobble_abond,]$aa_name,]$codon
  GC2_obs = rowSums(codon_usage[amino_acid_gc2[substr(amino_acid_gc2,3,3) %in% c("G","C")]],na.rm = T)
  ATGC2_obs = rowSums(codon_usage[amino_acid_gc2],na.rm = T)
  
  amino_acid_gc4 = tRNA_optimal[tRNA_optimal$nb_syn == 4,]$codon
  GC4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) %in% c("G","C")]],na.rm = T)
  ATGC4_obs = rowSums(codon_usage[amino_acid_gc4],na.rm = T)
  
  GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
  ATGC3_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)
  
  GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
  ATGCi_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)
  
  gc3 = (GC3_obs / ATGC3_obs)
  gci = (GCi_obs / ATGCi_obs)
  model_lm = lm( gc3 ~  gci)
  
  spearman_method_gc3gci = cor.test( gc3, gci,method = "spearman",exact=F)
  
  xaxis = codon_usage$median_fpkm 
  proportion = 5/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  ## Translational selection signal
  data_optiplus = data.frame()
  # for ( amino_acid in c("Ala", "Arg", "Gln", "Glu", "Gly", "Leu", "Lys", "Pro", "Ser", "Thr", "Val")){
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
  data_optiplus$rank = unlist(tapply(data_optiplus$expressed_overused_background,data_optiplus$amino_acid,function(x) rev(rank(x))))
  
  
  DUC_AG = data_optiplus[data_optiplus$expressed_overused_background > 0 ,]$codon
  DUC_AG = data_optiplus[data_optiplus$rank == 1 &
                           data_optiplus$codon %in% tRNA_optimal[tRNA_optimal$nb_syn >= 2,]$codon & 
                           !data_optiplus$amino_acid %in% c("Phe","Asn","Asp","His","Cys","Tyr","Met","Trp"),]$codon
  # DUC_AG = data_optiplus[!substr(data_optiplus$codon,3,3) %in% c("C","T"),]
  # DUC_AG = DUC_AG[order(DUC_AG$rank,decreasing = F),]
  # DUC_AG = DUC_AG[!duplicated(paste(DUC_AG$aa_name_scu)),]$codon
  
  
  
  dt_species = rbind(dt_species,data.frame(
    species,
    median_copy_tRNA = median(tRNA_optimal$nb_tRNA_copies / nb_genes),
    nb_genes,
    nb_genes_filtered = nrow(codon_usage),
    nb_codon_not_decoded,
    
    rho_aa_fpkm = spearman_method_aa$estimate,
    pval_aa_fpkm = spearman_method_aa$p.value,
    
    
    rho_gc3_gci = spearman_method_gc3gci$estimate,
    pval_gc3_gci = spearman_method_gc3gci$p.value,
    
    gc_duc_ag = sum(substr(DUC_AG,3,3) %in% c("G","C")) / length(DUC_AG),
    gc_abond_ag = sum(substr(abond_AG,3,3) %in% c("G","C")) / length(abond_AG),
    gc_abond_duet = sum(substr(abond_duet,3,3) %in% c("G","C")) / length(abond_duet),
    gc2 = mean(GC2_obs / ATGC2_obs,na.rm=T),
    gc2_top5 = tapply( GC2_obs / ATGC2_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    gc4 = mean(GC4_obs / ATGC4_obs,na.rm=T),
    gc4_top5 = tapply( GC4_obs / ATGC4_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    gc3 = mean(GC3_obs / ATGC3_obs,na.rm=T),
    gc3_top5 = tapply( GC3_obs / ATGC3_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    std_gc3 = stderror(GC3_obs / ATGC3_obs),
    var_gc3 = var(GC3_obs / ATGC3_obs,na.rm=T),
    
    gci = mean(GCi_obs / ATGCi_obs,na.rm=T),
    gci_top5 = tapply( GCi_obs / ATGCi_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    std_gci = stderror(GCi_obs / ATGCi_obs),
    var_gci = var(GCi_obs / ATGCi_obs,na.rm=T)
    
  ) )
  
  ### Translational selection intensity
  
  list_COA_global = c()
  list_COAneg_global = c()
  for ( type_aa in c( "WB_WC_notambiguous","WC_duet_ambiguous","global")){
    if (type_aa == "WB_WC_notambiguous" ){
      dt_selected = tRNA_optimal[ tRNA_optimal$nb_syn >= 2,]
      optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
      
      list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
      print(list_aa)
      list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon
      list_COA_neg = code[code$aa_name %in% list_aa,]$codon
      
      list_COA_global = append(list_COA_global,list_COA)
      list_COAneg_global = append(list_COAneg_global,list_COA_neg)
    } else if  (type_aa == "WC_duet_ambiguous" ){
      dt_selected = tRNA_optimal[  tRNA_optimal$aa_name  %in% c("Phe","Asn","Asp","His","Cys","Tyr") ,]
      optimal_count = table( dt_selected[dt_selected$Wobble_abond,]$aa_name )
      
      list_aa = names(optimal_count)[optimal_count != table(code$aa_name)[names(optimal_count)]]
      print(list_aa)
      list_COA = dt_selected[dt_selected$WC_abond & dt_selected$aa_name %in% list_aa,]$codon
      list_COA_neg = code[code$aa_name %in% list_aa ,]$codon
      list_COA_global = append(list_COA_global,list_COA)
      list_COAneg_global = append(list_COAneg_global,list_COA_neg)
    } else if  (type_aa == "global" ){
      
      list_COA = list_COA_global
      list_COA_neg = list_COAneg_global
    }
    print(list_COAneg_global)
    print(list_COA_global)
    
    
    if ( length(list_COA) != 0 ){
      
      ##### Over-used of pOC in expressed genes
      
      COA_obs = rowSums(codon_usage[ list_COA ],na.rm = T)
      COA_neg_obs = rowSums(codon_usage[ list_COA_neg ],na.rm = T)
      COA_obs_intronic = rowSums(codon_usage[paste(list_COA,'_intronic',sep = "")],na.rm = T)
      COA_neg_obs_intronic = rowSums(codon_usage[paste(list_COA_neg,'_intronic',sep = "")],na.rm = T)
      
      
      ## Over-used of pOC in constraint sites
      
      if (GTDrift_list_species[species,]$clade_group %in% c("Mammalia","Aves","Other Tetrapods")){
        data_conservation_sub = data_conservation_rmfirst1000bp[data_conservation_rmfirst1000bp$species == species,] 
      } else {
        data_conservation_sub = data_conservation[data_conservation$species == species,] 
      }
      
      table_constrain = data.frame(busco_id = data_conservation_sub$busco_id)
      for (constrain in c("_highconst","_modconst","_sligconst","_unconst")){
        table_constrain[,paste("COA",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_COA,constrain,sep = "")],na.rm = T)
        table_constrain[,paste("COA_neg",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_COA_neg,constrain,sep = "")],na.rm = T)
      }
      
      dt_translational_selection = data.frame(
        average_RTF_abundant_tRNA = mean( tRNA_optimal[list_COA,]$RTF[tRNA_optimal[list_COA,]$RTF != 0] ) ,
        expressed_overused = 100*(round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) - 
                                    round( mean( (COA_obs / COA_neg_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T) , 5)) ,
        expressed_overused_background = 100*((round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) - 
                                                round( mean( (COA_obs / COA_neg_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T) , 5)) - (
                                                  round(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) -
                                                    round( mean( (COA_obs_intronic / COA_neg_obs_intronic)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T),5)
                                                )),
        constraint_overused = 100*(mean(table_constrain$COA_highconst/table_constrain$COA_neg_highconst,na.rm = T) - 
                                     mean(table_constrain$COA_unconst/table_constrain$COA_neg_unconst,na.rm = T))
      )
    } else {
      dt_translational_selection = data.frame(
        average_RTF_abundant_tRNA = NA,
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
write.table(data1,"data/data1_bis.tab",quote=F,row.names = F,sep="\t")
