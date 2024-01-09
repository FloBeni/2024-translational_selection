# Generate Data 1

library(stringi)

path_data = "/home/fbenitiere/data/"

busco_tab = read.delim("/home/fbenitiere/data/Projet-SplicedVariants/DnDs/Metazoa_clades_v2/gene_No_aas_cds")
rownames(busco_tab) = busco_tab$species

data_conservation = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75.tab",sep=""))
data_conservation_rmfirst1000bp = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab",sep=""))

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

data1 = data.frame()
for (species in GTDrift_list_species$species){
  dt_species = data.frame()
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  if (
    file.exists(paste(path,"/tRNA_from_GFF.tab.gz",sep="")) &
    file.size(paste(path,"/tRNA_from_GFF.tab.gz",sep="")) != 38 ){
    tRNA_GFF = T
  } else if (
    file.exists(paste(path,"/tRNAscan_SE.tab.gz",sep="")) &
    file.size(paste(path,"/tRNAscan_SE.tab.gz",sep="")) != 36  ){
    tRNA_GFF = F
  } 
  
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
  abond_AG = tRNA_optimal[!tRNA_optimal$aa_name %in% c("Phe","Asn","Asp","His","Cys","Tyr","Met","Trp","Ile"),]
  abond_AG = abond_AG[substr(abond_AG$codon,3,3) %in% c("A","G"),]
  abond_AG = abond_AG[order(abond_AG$nb_tRNA_copies,decreasing = T),]
  abond_AG = abond_AG[ !abond_AG$aa_name_scu %in% abond_AG[abond_AG$nb_tRNA_copies == 0,]$aa_name_scu,]
  vector = (tapply(abond_AG$nb_tRNA_copies , abond_AG$aa_name_scu,function(x) sum(x == max(x))))
  abond_AG = abond_AG[abond_AG$aa_name_scu %in% names(vector[vector != 2]),]
  abond_AG = abond_AG[!duplicated(paste(abond_AG$aa_name_scu)),]$codon
  
  abundant_tRNA = tRNA_optimal[tRNA_optimal$nb_syn > 2,]
  abundant_tRNA = abundant_tRNA[abundant_tRNA$WC_abond,]$codon
  
  amino_acid_gc2 = tRNA_optimal[tRNA_optimal$nb_syn == 2 ,]$codon
  GC2_obs = rowSums(codon_usage[amino_acid_gc2[substr(amino_acid_gc2,3,3) %in% c("G","C")]],na.rm = T)
  ATGC2_obs = rowSums(codon_usage[amino_acid_gc2],na.rm = T)
  
  amino_acid_gc4 = tRNA_optimal[tRNA_optimal$nb_syn == 4,]$codon
  GC4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) %in% c("G","C")]],na.rm = T)
  ATGC4_obs = rowSums(codon_usage[amino_acid_gc4],na.rm = T)
  
  G4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) == "G"]] , na.rm = T)
  GC4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) %in% c("G","C")]] , na.rm = T)
  CT4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) %in% c("T","C")]] , na.rm = T)
  
  A4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) == "A"]] , na.rm = T)
  AT4_obs = rowSums(codon_usage[amino_acid_gc4[substr(amino_acid_gc4,3,3) %in% c("T","A")]] , na.rm = T)
  
  CTi_obs = rowSums(codon_usage[c("Ci","Ti")],na.rm = T)
  
  CT3_obs = rowSums(codon_usage[c("C3","T3")],na.rm = T)
  
  GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
  ATGCi_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)
  
  GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
  ATGC3_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)
  
  Gi_obs = rowSums(codon_usage[c("Gi")],na.rm = T)
  GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
  
  Ai_obs = rowSums(codon_usage[c("Ai")],na.rm = T)
  ATi_obs = rowSums(codon_usage[c("Ai","Ti")],na.rm = T)
  
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
  
  
  DUC = data_optiplus[data_optiplus$expressed_overused_background > 0 ,]$codon
  DUC = data_optiplus[data_optiplus$rank == 1 &
                        data_optiplus$codon %in% tRNA_optimal[tRNA_optimal$nb_syn >= 2,]$codon & 
                        !data_optiplus$amino_acid %in% c("Met","Trp"),]$codon
  
  
  ic_decoded = tRNA_optimal[substr(tRNA_optimal$codon,3,3) == "C" & tRNA_optimal$decoded & tRNA_optimal$nb_tRNA_copies == 0,]$codon
  
  DUC_IC = data_optiplus[data_optiplus$codon %in% c(ic_decoded , paste(sep="",substr(ic_decoded,1,2),"T")),]
  DUC_IC = DUC_IC[!duplicated(DUC_IC$amino_acid),]$codon
  
  
  gu_decoded = tRNA_optimal[substr(tRNA_optimal$codon,3,3) == "T" & tRNA_optimal$decoded & tRNA_optimal$nb_tRNA_copies == 0,]$codon
  
  DUC_GU = data_optiplus[data_optiplus$codon %in% c(gu_decoded , paste(sep="",substr(gu_decoded,1,2),"C")),]
  DUC_GU = DUC_GU[!duplicated(DUC_GU$amino_acid),]$codon
    
  dt_species = rbind(dt_species,data.frame(
    species,
    tRNA_GFF,
    prop_cds_expressed,
    # median_copy_tRNA = median(tRNA_optimal$nb_tRNA_copies / nb_genes),
    median_copy_tRNA_ponderate = median(tRNA_optimal$nb_tRNA_copies )/ nb_genes,
    median_copy_tRNA = median(tRNA_optimal$nb_tRNA_copies ),
    nb_genes,
    nb_genes_filtered = nrow(codon_usage),
    nb_codon_not_decoded,
    
    rho_aa_fpkm = spearman_method_aa$estimate,
    pval_aa_fpkm = spearman_method_aa$p.value,
    
    rho_gc3_gci = spearman_method_gc3gci$estimate,
    pval_gc3_gci = spearman_method_gc3gci$p.value,
    
    gc_duc = sum(substr(DUC,3,3) %in% c("G","C")) / length(DUC),
    g_abond_ag = sum(substr(abond_AG,3,3) %in% c("G")) / length(abond_AG),
    c_duc_ic = sum(substr(DUC_IC,3,3) %in% c("C")) / length(DUC_IC),
    c_duc_gu = sum(substr(DUC_GU,3,3) %in% c("C")) / length(DUC_GU),
    ct_abundant_tRNA4 = sum(substr(abundant_tRNA,3,3) %in% c("C","T")) / length(abundant_tRNA),
    
    gc3_expression = sum(observation[c("C3","G3")],na.rm = T) / sum(observation[c("A3","T3","C3","G3")],na.rm = T),
    ct3 = mean(CT3_obs / ATGC3_obs,na.rm=T),
    gc3 = mean(GC3_obs / ATGC3_obs,na.rm=T),
    gc3_top5 = tapply( GC3_obs / ATGC3_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    std_gc3 = stderror(GC3_obs / ATGC3_obs),
    var_gc3 = var(GC3_obs / ATGC3_obs,na.rm=T),
    
    g4 = mean(G4_obs / GC4_obs,na.rm=T),
    g4_top5 = tapply( G4_obs / GC4_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    a4 = mean(A4_obs / AT4_obs,na.rm=T),
    a4_top5 = tapply( A4_obs / AT4_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    ct4 = mean(CT4_obs / ATGC4_obs,na.rm=T),
    gc4 = mean(GC4_obs / ATGC4_obs,na.rm=T),
    gc4_top5 = tapply( GC4_obs / ATGC4_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    gc2 = mean(GC2_obs / ATGC2_obs,na.rm=T),
    gc2_top5 = tapply( GC2_obs / ATGC2_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    
    gi = mean(Gi_obs / GCi_obs,na.rm=T),
    gi_top5 = tapply( Gi_obs / GCi_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    ai = mean(Ai_obs / ATi_obs,na.rm=T),
    ai_top5 = tapply( Ai_obs / ATi_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
    cti = mean(CTi_obs / ATGCi_obs,na.rm=T),
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

      ##### Over-used of POC in expressed genes

      COA_obs = rowSums(codon_usage[ list_COA ],na.rm = T)
      COA_neg_obs = rowSums(codon_usage[ list_COA_neg ],na.rm = T)
      COA_obs_intronic = rowSums(codon_usage[paste(list_COA,'_intronic',sep = "")],na.rm = T)
      COA_neg_obs_intronic = rowSums(codon_usage[paste(list_COA_neg,'_intronic',sep = "")],na.rm = T)


      ## Over-used of POC in constraint sites

      if (GTDrift_list_species[species,]$clade_group %in% c("Mammalia","Aves","Other Tetrapods")){
        data_conservation_sub = data_conservation_rmfirst1000bp[data_conservation_rmfirst1000bp$species == species & data_conservation_rmfirst1000bp$protein %in% codon_usage$protein_id,]
      } else {
        data_conservation_sub = data_conservation[data_conservation$species == species ,]
      }


      table_constrain = data.frame(busco_id = data_conservation_sub$busco_id)
      for (constrain in c("_highconst","_modconst","_sligconst","_unconst")){
        table_constrain[,paste("COA",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_COA,constrain,sep = "")],na.rm = T)
        table_constrain[,paste("COA_neg",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_COA_neg,constrain,sep = "")],na.rm = T)
      }

      dt_translational_selection = data.frame(
        nb_aa = length(list_aa),
        nb_busco = busco_tab[species,]$No_CDS_Busco,
        average_RTF_abundant_tRNA = mean( tRNA_optimal[list_COA,]$RTF[tRNA_optimal[list_COA,]$WC_abond] ) ,
        average_RTFprime_abundant_tRNA = mean( tRNA_optimal[list_COA,]$RTF_prime[tRNA_optimal[list_COA,]$WC_abond ] ) ,
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
        nb_aa = NA,
        nb_busco = NA,
        average_RTF_abundant_tRNA = NA,
        average_RTFprime_abundant_tRNA = NA,
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
