# Generate Data 13
path = "/home/fbenitiere/data/"

code = read.delim(paste(path,"Projet-SplicedVariants/Fichiers-data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )




amino_acid_list = unique(code[code$nb_syn > 1 & code$aa_name != "Ter",]$aa_name)


clade_dt = read.delim(paste( "data/clade_dt.tab",sep=""),header=T)
rownames(clade_dt) = clade_dt$species
clade_dt$clade_group = factor(clade_dt$clade_group, levels = c("Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Teleostei","Mammalia","Aves","Other Tetrapodes"))

clade_dt = clade_dt[ clade_dt$clade_group %in% c("Lepido Diptera") & clade_dt$species %in% list_species,]


species = "Drosophila_melanogaster"
species = "Eumeta_japonica"
species = "Lucilia_cuprina"
species = "Melitaea_cinxia"

data13 = data.frame()
for (species in clade_dt$species){print(species)
  
  codon_usage = read.delim( paste(path,"/Projet-SplicedVariants/Analyses/",species,"/codon_usage_gene_fpkm.tab",sep="") )
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  codon_usage = codon_usage[codon_usage$length_cds == codon_usage$length,]
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage$stop_codon = rowSums(codon_usage[,stop_codon])
  if (!species %in% c("Melospiza_melodia","Oenanthe_oenanthe","Eudyptes_filholi")){codon_usage = codon_usage[codon_usage$stop_codon == 1,]}
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 & !is.na(codon_usage$median_fpkm) ,]
  
  
  
  code_table = read.delim(paste(path,"/Projet-NeGA/translational_selection/ExpOpti/" , species , "_code_table.tab" , sep=""))
  rownames(code_table) = code_table$codon
  code_table = code_table[code_table$aa_name %in% amino_acid_list,]
  
  xaxis = codon_usage$median_fpkm 
  proportion = 5/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  FPKM_bins = log10(FPKM_bins)
  
  
  data_optiplus = data.frame()
  
  for ( amino_acid in amino_acid_list){
    codon_used = rownames(code[code$aa_name == amino_acid,])
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
        optifreq_top5 = round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5),
        opti_freq_low50 = round( mean( COA_obs[codon_usage$median_fpkm < median(codon_usage$median_fpkm )] /
                                         COA_neg_obs[codon_usage$median_fpkm < median(codon_usage$median_fpkm )] , na.rm=T) , 5),
        
        optifreq_top5_intron = round(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5),
        opti_freq_low50_intron = round( mean( COA_obs_intronic[codon_usage$median_fpkm < median(codon_usage$median_fpkm )] /
                                                COA_neg_obs_intronic[codon_usage$median_fpkm < median(codon_usage$median_fpkm )] , na.rm=T),5)
        
      ))
    }
  }
  
  data_optiplus$ecart = (data_optiplus$optifreq_top5 - data_optiplus$optifreq_top5_intron) - (data_optiplus$opti_freq_low50 - data_optiplus$opti_freq_low50_intron)
  
  data_optiplus = data_optiplus[order(data_optiplus$ecart,decreasing = T),]
  
  code_table = read.delim(paste(path,"/Projet-NeGA/translational_selection/ExpOpti/" , species , "_code_table.tab" , sep=""))
  dt_selected = code_table[ code_table$nb_syn >= 2,]
  optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )

  list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
  nb_aa = length(list_aa)

  list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon

  data_optiplus$list_COA = data_optiplus$codon %in% list_COA 
  
  data_optiplus = data_optiplus[ !duplicated( data_optiplus$amino_acid ) , ]
  print( table( data_optiplus$ecart > 0 ) )
  # data_optiplus = data_optiplus[ data_optiplus$ecart > 0 , ]
  print(table(data_optiplus[data_optiplus$amino_acid %in%list_aa,]$codon %in% list_COA))
  data_optiplus_quart = data_optiplus[data_optiplus$amino_acid %in% code_table[code_table$nb_syn == 4 ,]$aa_name,]
  
  data13 = rbind(data13,data.frame(species,
                                   mean_ecart = mean(data_optiplus$ecart),
                                   prop_abond = sum(data_optiplus[data_optiplus$amino_acid %in% list_aa,]$codon %in% list_COA) / length( data_optiplus[data_optiplus$amino_acid %in% list_aa,]$codon ),
                                   GC_DUC =  sum(substr(data_optiplus$codon , 3 , 3) %in% c("G" , "C")) / length(data_optiplus$codon),
                                   mean_ecart_quart = mean(data_optiplus_quart$ecart),
                                   prop_abond_quart = sum(data_optiplus_quart[data_optiplus_quart$amino_acid %in% list_aa,]$codon %in% list_COA) / length( data_optiplus_quart[data_optiplus_quart$amino_acid %in%list_aa,]$codon ),
                                   GC_DUC_quart =  sum(substr(data_optiplus_quart$codon , 3 , 3) %in% c("G" , "C")) / length(data_optiplus_quart$codon)
  ))
  
}


write.table(data13 , "data/data13.tab" , quote=F , row.names = F , sep="\t")


