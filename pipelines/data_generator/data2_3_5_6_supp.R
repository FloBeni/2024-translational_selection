# Generate Data 2, 3, 5 and 6

code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])

data_conservation = read.delim(paste("data/compilation_prop_gap_pergene_25_50_75.tab.gz",sep=""),comment.char = "#")
data_conservation_rmfirst1000bp = read.delim(paste("data/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab.gz",sep=""),comment.char = "#")

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab",comment.char = "#")
rownames(GTDrift_list_species) = GTDrift_list_species$species


species_list = c( "Caenorhabditis_elegans" , "Drosophila_melanogaster" , "Homo_sapiens" , "Musca_domestica" ,"Anopheles_gambiae","Pieris_rapae",
                  "Bactrocera_oleae","Ceratitis_capitata","Bombus_terrestris","Ignelater_luminosus" )

data2 = data.frame()
data3 = data.frame()
data5 = data.frame()
data6 = data.frame()
for (species in species_list){
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.txt.gz",sep=""),comment.char = "#")
  
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
    gene_id = codon_usage$gene_id,
    GCi,
    GC3
  ))
  
  
  ##### AA frequency
  
  observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )
  
  
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""),comment.char = "#")
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
                      tRNA_gene_copy= sum(tRNA_optimal[codon_used,]$nb_tRNA_copies,na.rm = T),
                      obs_aminoacid = sum(observation[codon_used])
                    ))
  }
  aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]
  
  
  
  data3 = rbind(data3,aa_data)
  
  
  ##### Over-used of pOC in expressed genes
  
  xaxis = codon_usage$median_fpkm 
  proportion = 2/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  for ( type_aa in c( "POC1","POC2","POCs")){
    if (type_aa == "POC1" ){
      subset_selected = tRNA_optimal[tRNA_optimal$POC1,]
      list_POC = subset_selected$codon
      list_aa = subset_selected$aa_name
      list_codon = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa,]$codon
    } else if  (type_aa == "POC2" ){
      subset_selected = tRNA_optimal[tRNA_optimal$POC2,]
      list_POC = subset_selected$codon
      list_aa = subset_selected$aa_name
      list_codon = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa,]$codon
    } else if  (type_aa == "POCs" ){
      subset_selected = tRNA_optimal[tRNA_optimal$POC2 | tRNA_optimal$POC1,]
      list_POC = subset_selected$codon
      list_aa = subset_selected$aa_name
      list_codon = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa,]$codon
    }
    print(list_POC)
    print(list_codon)
    
    
    
    POC_obs = rowSums(codon_usage[ list_POC ],na.rm = T)
    POC_codon = rowSums(codon_usage[ list_codon ],na.rm = T)
    
    POC_obs_intronic = rowSums(codon_usage[paste(list_POC,'_intronic',sep = "")],na.rm = T)
    POC_codon_intronic = rowSums(codon_usage[paste(list_codon,'_intronic',sep = "")],na.rm = T)
    
    data5 = rbind(data5,data.frame(
      species,
      set = type_aa,
      nb_genes = as.numeric(table(intervalle_FPKM)),
      freq =  tapply( POC_obs / POC_codon  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
      fpkm = FPKM_bins,
      categorie = "Putative optimal codons (POC)",
      type_aa=paste(type_aa," codons, Nb aa = ",length(list_aa),sep=""),
      gene_set = paste("N = " , format(sum(!is.na(POC_obs / POC_codon) ),big.mark=",",scientific=T),sep=""),
      nb_poc =  length(list_POC) ,
      nb_aa = length(list_aa),
      list=paste(list_POC,collapse =";")
    ))
    
    
    data5 = rbind(data5,data.frame(
      species,
      set = type_aa,
      nb_genes = as.numeric(table(intervalle_FPKM)),
      freq =   tapply( POC_obs_intronic / POC_codon_intronic  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
      fpkm = FPKM_bins,
      categorie = "POC-matching triplets (POCMT)",
      type_aa=paste(type_aa," codons, Nb aa = ",length(list_aa),sep=""),
      gene_set=paste("N = " , format(sum(!is.na(POC_obs_intronic / POC_codon_intronic) ),big.mark=",",scientific=T),sep=""),
      nb_poc =  length(list_POC),
      nb_aa = length(list_aa),
      list=paste(list_POC,collapse=";")
    ))
    
    
    ## Over-used of POC in constraint sites
    
    if (GTDrift_list_species[species,]$clade_group %in% c("Mammalia","Aves","Other Tetrapods")){
      data_conservation_sub = data_conservation_rmfirst1000bp[data_conservation_rmfirst1000bp$species == species & data_conservation_rmfirst1000bp$protein %in% codon_usage$protein_id,] 
    } else {
      data_conservation_sub = data_conservation[data_conservation$species == species & data_conservation$protein %in% codon_usage$protein_id,] 
    }
    
    table_constrain = data.frame(busco_id = data_conservation_sub$busco_id)
    for (constrain in c("_highconst","_modconst","_sligconst","_unconst")){
      table_constrain[,paste("POC",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_POC,constrain,sep = "")],na.rm = T)
      table_constrain[,paste("POC_codon",constrain,sep="")] = rowSums(data_conservation_sub[paste(list_codon,constrain,sep = "")],na.rm = T)
    }
    
    data6 = rbind(data6, data.frame(
      species,
      categorie = "Highly constrained",
      set = type_aa,
      busco_id=table_constrain$busco_id,
      freq = table_constrain$POC_highconst / table_constrain$POC_codon_highconst ,
      nb_site = sum(data_conservation_sub$len_high_const_seq),
      nb_genes = nrow(data_conservation_sub)
    ))
    
    data6 = rbind(data6, data.frame(
      species,
      categorie = "Moderately constrained",
      set = type_aa,
      busco_id=table_constrain$busco_id,
      freq = table_constrain$POC_modconst / table_constrain$POC_codon_modconst ,
      nb_site = sum(data_conservation_sub$len_mod_const_seq),
      nb_genes = nrow(data_conservation_sub)
    ))
    data6 = rbind(data6, data.frame(
      species,
      categorie = "Slighlty constrained",
      set = type_aa,
      busco_id=table_constrain$busco_id,
      freq = table_constrain$POC_sligconst / table_constrain$POC_codon_sligconst ,
      nb_site = sum(data_conservation_sub$len_slight_const_seq),
      nb_genes = nrow(data_conservation_sub)
    ))
    data6 = rbind(data6, data.frame(
      species,
      categorie = "Unconstrained",
      set = type_aa,
      busco_id=table_constrain$busco_id,
      freq = table_constrain$POC_unconst / table_constrain$POC_codon_unconst ,
      nb_site = sum(data_conservation_sub$len_unconst_seq),
      nb_genes = nrow(data_conservation_sub)
    ))
  }
}

write.table(data2 , "data/data2_supp.tab",quote=F,row.names = F,sep="\t")
write.table(data3 , "data/data3_supp.tab",quote=F,row.names = F,sep="\t")
write.table(data5 , "data/data5_supp.tab",quote=F,row.names = F,sep="\t")
write.table(data6 , "data/data6_supp.tab",quote=F,row.names = F,sep="\t")
