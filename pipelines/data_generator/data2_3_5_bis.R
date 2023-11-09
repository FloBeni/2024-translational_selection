# Generate Data 2, 3 and 5

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


species_list = c("Caenorhabditis_elegans","Drosophila_melanogaster","Homo_sapiens","Musca_domestica")

data2 = data.frame()
data3 = data.frame()
data5 = data.frame()
for (species in species_list){
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxID",taxID,"/",genome_assembly,sep="")
  
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
  
  if (length(list_COA) != 0){
    COA_obs_intronic = rowSums(codon_usage[paste(list_COA,'_intronic',sep = "")],na.rm = T)
    COA_neg_obs_intronic = rowSums(codon_usage[paste(list_COA_neg,'_intronic',sep = "")],na.rm = T)
  } else {
    COA_obs_intronic = rowSums(codon_usage_selected[list_COA],na.rm = T)
    COA_neg_obs_intronic = rowSums(codon_usage_selected[list_COA_neg],na.rm = T)
  }
  
  
  data5 = rbind(data5,data.frame(
    species,
    freq =  tapply( COA_obs / COA_neg_obs  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
    fpkm = FPKM_bins,
    categorie = "optimal codons",
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
    categorie = "intronic triplets control",
    type_aa=paste("Wb_WC_notambiguous"," codons, Nb aa = ",nb_aa,sep=""),
    set=paste("N = " , sum(ATGCi_neg_obs != 0),sep=""),
    nb_poc ,
    nb_aa,
    list=paste(list_COA,collapse=";")
  ))
  
}

write.table(data2 , "data/data2_bis.tab",quote=F,row.names = F,sep="\t")
write.table(data3 , "data/data3_bis.tab",quote=F,row.names = F,sep="\t")
write.table(data5 , "data/data5_bis.tab",quote=F,row.names = F,sep="\t")
