# Generate Data 4
library(stringi)


wobble_pairing = c("C"="IC","T"="GU","G"="UG","A"="IA")

wobble_associat_wc = list("T"="C",
                          "A"="A",
                          "G"="A",
                          "C"="T")


code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")
stop_codon = rownames(code[code$aa_name == "Ter",])


amino_acid_list = unique(code[code$nb_syn > 1 & code$aa_name != "Ter",]$aa_name_scu)
amino_acid_list = unique(code[ code$aa_name != "Ter",]$aa_name)

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


list_species = c("Homo_sapiens","Caenorhabditis_elegans","Drosophila_melanogaster","Musca_domestica")


data4 = data.frame()
for (species in list_species){
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
  
  
  if( nrow(codon_usage) == 0){filtering = "not enough CDS"  } else {
    if (
      file.exists(paste(path,"/tRNAscan.tab.gz",sep="")) &
      file.size(paste(path,"/tRNAscan.tab.gz",sep="")) != 33 ){
      tRNASE_gff = read.delim( paste(path,"/tRNAscan.tab.gz",sep="") )
      tRNASE_gff$codon = sapply(tRNASE_gff$anticodon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
      tRNASE_gff_table = table(tRNASE_gff$codon)
      tRNASE_copies_table = tRNASE_gff_table
      tRNA_GFF = "from GFF"
    } else if (
      file.exists(paste(path,"/tRNAscan_SE.tab.gz",sep="")) &
      file.size(paste(path,"/tRNAscan_SE.tab.gz",sep="")) != 33  ){
      tRNASE_copies = read.delim(paste(path,"/tRNAscan_SE.tab.gz",sep=""), header = F)
      colnames(tRNASE_copies) = unlist(tRNASE_copies[2,])
      tRNASE_copies = tRNASE_copies[4:nrow(tRNASE_copies),]
      tRNASE_copies = tRNASE_copies[as.vector(tRNASE_copies$Note) != "pseudo",]
      tRNASE_copies = tRNASE_copies[tRNASE_copies$Note != "pseudo",]
      tRNASE_copies = tRNASE_copies[as.numeric(tRNASE_copies$Score) > 55 & !is.na(as.numeric(tRNASE_copies$Score)),]
      tRNASE_copies$anticodon = sapply(tRNASE_copies$Codon, function(x) chartr("TUACG","TTACG",x)  )
      tRNASE_copies$codon = sapply(tRNASE_copies$Codon, function(x) chartr("TUACG","AATGC",stri_reverse(x)) )
      tRNASE_copies_table = table(tRNASE_copies$codon)
      tRNA_GFF = "not from GFF"
    } else { tRNASE_copies_table = 0
    tRNA_GFF = "" }
    
    
    codon_usage$gene_length = codon_usage$length_cds
    
    code_table = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
    code_table$wb = wobble_pairing[substr(code_table$codon,3,3)]
    code_table$decoded_codon = sapply(code_table$codon,function(x) paste(substr(x,1,2),wobble_associat_wc[[substr(x,3,3)]],collapse="",sep=""))
    rownames(code_table) = code_table$codon
    code_table = code_table[code_table$aa_name %in% amino_acid_list,]
    nb_codon_not_decoded = sum(!code_table$decoded)
    
    
    nb_long_length = sum(codon_usage$length_cds >= median(codon_usage$length_cds))
    nb_short_length = sum(codon_usage$length_cds < median(codon_usage$length_cds))
    
    vector_count = c("long"= nb_long_length , "short"=nb_short_length, "all"=nb_short_length+nb_long_length)
    
    length_limit = "all"
    codon_usage_selected = codon_usage
    
    
    xaxis = codon_usage_selected$median_fpkm 
    proportion = 5/100
    quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
    intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
    
    FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
    prop_GU_or_IC = NA
    # for ( type_aa in c("IC","GU",amino_acid_list,c("Wb_WC_notambiguous","WC_duet","WC_duet_ambiguous","WC_duet_notambiguous"))){
    for ( type_aa in c("Wb_WC_notambiguous","WC_duet_ambiguous","WB_vs_WC_notduet")){
      if (type_aa == "WC_duet_notambiguous" ){
        dt_selected = code_table[ code_table$nb_syn == 2,]
        optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
        
        list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[ dt_selected$WC_abond & dt_selected$aa_name %in% list_aa,]$codon
        list_COA_neg = code[code$aa_name %in% list_aa,]$codon
        nb_codon = length(list_COA)
        
        
      }  else if (type_aa == "Wb_WC_notambiguous" ){
        dt_selected = code_table[ code_table$nb_syn >= 2,]
        optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
        
        list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon
        list_COA_neg = code[code$aa_name %in% list_aa,]$codon
        nb_codon = length(list_COA)
        
        
      }  else if (type_aa == "WC_duet" ){
        dt_selected = code_table[ code_table$nb_syn == 2,]
        optimal_count = table( dt_selected[ dt_selected$WC_abond,]$aa_name )
        
        list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[ dt_selected$WC_abond & dt_selected$aa_name %in% list_aa,]$codon
        list_COA_neg = code[code$aa_name %in% list_aa,]$codon
        nb_codon = length(list_COA)
        
      }  else if  (type_aa == "WC_duet_ambiguous" ){
        dt_selected = code_table[  code_table$nb_syn  == 2 ,]
        optimal_count = table( dt_selected[dt_selected$Wobble_abond,]$aa_name )
        
        list_aa = names(optimal_count)[optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[dt_selected$WC_abond & dt_selected$aa_name %in% list_aa,]$codon
        list_COA_neg = code[code$aa_name %in% list_aa ,]$codon
        nb_codon = length(list_COA)
        
        prop_GU_or_IC = sum(wobble_pairing[substr(dt_selected[dt_selected$Wobble_abond & dt_selected$aa_name %in% list_aa,]$codon,3,3)] == "GU")
        
      } else if  (type_aa %in% c("IC","GU" )){
        dt_selected = code_table[ code_table$nb_syn >= 2,]
        optimal_count = table( dt_selected[ dt_selected$Wobble_abond & dt_selected$wb == type_aa,]$aa_name )
        
        list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[ dt_selected$Wobble_abond & dt_selected$wb == type_aa & dt_selected$aa_name %in% list_aa,]$decoded_codon
        list_COA_neg =  c(dt_selected[ dt_selected$Wobble_abond & dt_selected$wb == type_aa & dt_selected$aa_name %in% list_aa,]$codon,list_COA)
        nb_codon = length(list_COA)
        
        
      } else if  (type_aa == "WB_vs_WC_notduet"){
        dt_selected = code_table[ code_table$nb_syn > 2,]
        optimal_count = table( dt_selected[ dt_selected$Wobble_abond ,]$aa_name )
        
        list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[ dt_selected$Wobble_abond & dt_selected$aa_name %in% list_aa,]$decoded_codon
        list_COA_neg =  c(dt_selected[ dt_selected$Wobble_abond & dt_selected$aa_name %in% list_aa,]$codon,list_COA)
        nb_codon = length(list_COA)
        print(wobble_pairing[substr(dt_selected[dt_selected$Wobble_abond & dt_selected$aa_name %in% list_aa,]$codon,3,3)])
        prop_GU_or_IC = sum( wobble_pairing[substr(dt_selected[dt_selected$Wobble_abond & dt_selected$aa_name %in% list_aa,]$codon,3,3)] == "IC")
        
      }  else {
        dt_selected = code_table[ code_table$aa_name == type_aa,]
        optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
        
        list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
        nb_aa = length(list_aa)
        
        list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon
        list_COA_neg = code[code$aa_name %in% list_aa,]$codon
        nb_codon = length(list_COA)
        
      }
      
      COA_obs = rowSums(codon_usage_selected[ list_COA ],na.rm = T)
      COA_neg_obs = rowSums(codon_usage_selected[ list_COA_neg ],na.rm = T)
      
      
      if (length(list_COA) != 0){
        COA_obs_intronic = rowSums(codon_usage_selected[paste(list_COA,'_intronic',sep = "")],na.rm = T)
        COA_neg_obs_intronic = rowSums(codon_usage_selected[paste(list_COA_neg,'_intronic',sep = "")],na.rm = T)
      } else {
        COA_obs_intronic = rowSums(codon_usage_selected[list_COA],na.rm = T)
        COA_neg_obs_intronic = rowSums(codon_usage_selected[list_COA_neg],na.rm = T)
      }
      
      
      data4 = rbind(data4,data.frame(
        species,
        freq =  tapply( COA_obs / COA_neg_obs  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
        fpkm = FPKM_bins,
        categorie = "optimal codons",
        group = "not_control",
        type_aa=paste(type_aa," codons, Nb aa = ",nb_aa,sep=""),
        set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
        nb_codon ,
        nb_aa,
        prop_GU_or_IC,
        list=paste(list_COA,collapse =";")
      ))
      
      
      data4 = rbind(data4,data.frame(
        species,
        freq =   tapply( COA_obs_intronic / COA_neg_obs_intronic  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
        fpkm = FPKM_bins,
        categorie = "intronic triplets control",
        group = "control",
        type_aa=paste(type_aa," codons, Nb aa = ",nb_aa,sep=""),
        set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
        nb_codon ,
        nb_aa,
        prop_GU_or_IC,
        list=paste(list_COA,collapse=";")
      ))
    }
  }
}

write.table(data4,"data/data4.tab",quote=F,row.names = F,sep="\t")

