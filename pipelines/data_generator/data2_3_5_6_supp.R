# Generate Data 2, 3, 5 and 6

path_data = "/home/fbenitiere/data/"

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])

data_conservation = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75.tab",sep=""))
data_conservation_rmfirst1000bp = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab",sep=""))

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

species_list = c( "Caenorhabditis_elegans" , "Drosophila_melanogaster" , "Homo_sapiens" , "Musca_domestica" )

data2 = data.frame()
data3 = data.frame()
data5 = data.frame()
data6 = data.frame()
for (species in species_list){
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )
  
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  
  codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]
  
  if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
    codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}
  
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
  
  ## GCi vs GC3
  
  GC3_obs = rowSums(codon_usage[c("C3","G3")],na.rm = T)
  ATGC3_neg_obs = rowSums(codon_usage[c("A3","T3","C3","G3")],na.rm = T)
  
  GCi_obs = rowSums(codon_usage[c("Ci","Gi")],na.rm = T)
  ATGCi_neg_obs = rowSums(codon_usage[c("Ai","Ti","Ci","Gi")],na.rm = T)
  
  GC3 = GC3_obs / ATGC3_neg_obs
  GCi = GCi_obs / ATGCi_neg_obs
  
  data2 = rbind(data2,data.frame(
    species,
    GCi,
    GC3
  ))
  
  
  ##### AA frequency
  
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
                      species,
                      amino_acid ,
                      letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                      tRNASE_copies= sum(tRNA_optimal[codon_used,]$nb_tRNA_copies,na.rm = T),
                      obs_codon = sum(observation[codon_used])
                    ))
  }
  aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]
  
  
  
  data3 = rbind(data3,aa_data)
  
  
  ##### Over-used of pOC in expressed genes
  
  xaxis = codon_usage$median_fpkm 
  proportion = 5/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  dt_selected = tRNA_optimal[ tRNA_optimal$nb_syn >= 2,]
  optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
  
  list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
  nb_aa = length(list_aa)
  
  list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon
  list_COA_neg = code[code$aa_name %in% list_aa,]$codon
  nb_poc = length(list_COA)
  
  COA_obs = rowSums(codon_usage[ list_COA ],na.rm = T)
  COA_neg_obs = rowSums(codon_usage[ list_COA_neg ],na.rm = T)
  
  COA_obs_intronic = rowSums(codon_usage[paste(list_COA,'_intronic',sep = "")],na.rm = T)
  COA_neg_obs_intronic = rowSums(codon_usage[paste(list_COA_neg,'_intronic',sep = "")],na.rm = T)
  
  data5 = rbind(data5,data.frame(
    species,
    freq =  tapply( COA_obs / COA_neg_obs  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
    fpkm = FPKM_bins,
    categorie = "Putative optimal codons (POC)",
    type_aa=paste("Wb_WC_notambiguous"," codons, Nb aa = ",nb_aa,sep=""),
    set=paste("N = " , sum(ATGC3_neg_obs != 0),sep=""),
    nb_poc ,
    nb_aa,
    list=paste(list_COA,collapse =";")
  ))
  
  
  data5 = rbind(data5,data.frame(
    species,
    freq =   tapply( COA_obs_intronic / COA_neg_obs_intronic  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
    fpkm = FPKM_bins,
    categorie = "POC-matching triplets (POCMT)",
    type_aa=paste("Wb_WC_notambiguous"," codons, Nb aa = ",nb_aa,sep=""),
    set=paste("N = " , sum(ATGCi_neg_obs != 0),sep=""),
    nb_poc ,
    nb_aa,
    list=paste(list_COA,collapse=";")
  ))
  
  
  ## Over-used of POC in constraint sites
  
  if (GTDrift_list_species[species,]$clade_group %in% c("Mammalia","Aves","Other Tetrapods")){
    data_conservation_sub = data_conservation_rmfirst1000bp[data_conservation_rmfirst1000bp$species == species & data_conservation_rmfirst1000bp$protein %in% codon_usage$protein_id,] 
  } else {
    data_conservation_sub = data_conservation[data_conservation$species == species & data_conservation$protein %in% codon_usage$protein_id,] 
  }
  
  table_constrain = data.frame(busco_id = data_conservation_sub$busco_id)
  for (constrain in c("_highconst","_modconst","_sligconst","_unconst")){
    table_constrain[,paste("COA",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_COA,constrain,sep = "")],na.rm = T)
    table_constrain[,paste("COA_neg",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_COA_neg,constrain,sep = "")],na.rm = T)
  }
  
  data6 = rbind(data6, data.frame(
    species,
    freq = table_constrain$COA_highconst / table_constrain$COA_neg_highconst ,
    nb_site = sum(data_conservation_sub$len_high_const_seq),
    nb_genes = nrow(data_conservation_sub),
    categorie = "Highly constrained"
  ))
  
  data6 = rbind(data6, data.frame(
    species,
    freq = table_constrain$COA_modconst / table_constrain$COA_neg_modconst ,
    nb_site = sum(data_conservation_sub$len_mod_const_seq),
    nb_genes = nrow(data_conservation_sub),
    categorie = "Moderately constrained"
  ))
  data6 = rbind(data6, data.frame(
    species,
    freq = table_constrain$COA_sligconst / table_constrain$COA_neg_sligconst ,
    nb_site = sum(data_conservation_sub$len_slight_const_seq),
    nb_genes = nrow(data_conservation_sub),
    categorie = "Slighlty constrained"
  ))
  data6 = rbind(data6, data.frame(
    species,
    freq = table_constrain$COA_unconst / table_constrain$COA_neg_unconst ,
    nb_site = sum(data_conservation_sub$len_unconst_seq),
    nb_genes = nrow(data_conservation_sub),
    categorie = "Unconstrained"
  ))
}

write.table(data2 , "data/data2_supp.tab",quote=F,row.names = F,sep="\t")
write.table(data3 , "data/data3_supp.tab",quote=F,row.names = F,sep="\t")
write.table(data5 , "data/data5_supp.tab",quote=F,row.names = F,sep="\t")
write.table(data6 , "data/data6_supp.tab",quote=F,row.names = F,sep="\t")
