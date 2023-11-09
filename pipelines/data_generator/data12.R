# Generate Data 12
path = "/home/fbenitiere/data/"
# path = "/beegfs/data/fbenitiere/"

library(stringi)
stderror <- function(x) sd(x , na.rm = T)/sqrt(length(x[!is.na(x)] ))


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


wobble_rule = list("T"=c("A","G"),
                   "A"=c("T","A"),
                   "G"=c("C","T"),
                   "C"=c("G","A"))

wobble_pairing = c("C"="IC","T"="GU","G"="UG","A"="IA")

wobble_associat_wc = list("T"="C",
                          "A"="A",
                          "G"="A",
                          "C"="T")



GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


data1 = read.delim("data/data1.tab")
data1 = data1[data1$pval_aa_fpkm < 0.05,]
data1 = data1[data1$nb_gene > 5000 & data1$pval_aa_fpkm < 0.05,]
# data1 = data1[data1$nb_gene >= 10000 ,]
list_species = unique(data1$species)

i = 0
data12 = data.frame()
for( species in list_species ){ print(species)
  i=i+1
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
  print(round(table(codon_usage$median_fpkm == 0)['TRUE']/sum(table(codon_usage$median_fpkm == 0))*100))
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
    
    if ( length(tRNASE_copies_table) == 1 ){ filtering = "no tRNA copies file"
    print("####### problem")
    
    } else {
      
      observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )
      # observation = colSums( codon_usage[3:70] ,na.rm = T )
      
      aa_data = data.frame()
      for (amino_acid in unique(code$aa_name)){
        codon_used = rownames(code[code$aa_name == amino_acid,])
        aa_data = rbind(aa_data,
                        data.frame(
                          amino_acid ,
                          letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                          tRNASE_copies= sum(tRNASE_copies_table[codon_used],na.rm = T),
                          obs_codon = sum(observation[codon_used])
                        ))
      }
      
      aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]
      
      spearman_method_aa = cor.test( aa_data$tRNASE_copies, aa_data$obs_codon,method="spearman",exact=F)
      
      ######## CODON OPTIMAL ########
      
      aa_data = data.frame()
      for (amino_acid in unique(code$aa_name)){
        codon_used = rownames(code[code$aa_name == amino_acid,])
        
        nb_syn = length(codon_used)
        nb_aa_copie = sum(tRNASE_copies_table[codon_used],na.rm = T)/nb_syn
        nb_observation = sum(observation[codon_used],na.rm = T)/nb_syn
        
        for (codon in codon_used){
          if( sum(tRNASE_copies_table[codon] , na.rm = T) == 0 ){
            anticodon = paste(unlist(wobble_rule[substr(codon,3,3)]),
                              chartr( "TUACG" , "AATGC" , stri_reverse(substr(codon,1,2))),sep = "")
            
            if ( length( code[ code$anticodon %in% anticodon & code$aa_name == amino_acid, ]$anticodon ) != length(anticodon) ){
              # print("GIGA PROBLEM")
              # print(amino_acid)
              # print(codon)
            }
            
            anticodon = code[  code$anticodon %in% anticodon & code$aa_name == amino_acid, ]$anticodon
            
            tRNASE_copies_table_antic = tRNASE_copies_table
            names(tRNASE_copies_table_antic) = sapply(names(tRNASE_copies_table_antic),function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
            
            aa_data = rbind(aa_data,
                            data.frame(
                              codon ,
                              nb_syn,
                              amino_acid,
                              rules = "wobble",
                              letter_aa = codon ,
                              anticodon = chartr("TUACG","AATGC",stri_reverse(codon)),
                              count = sum(observation[codon]) / sum(observation[codon_used],na.rm = T),
                              RTF = sum(tRNASE_copies_table_antic[anticodon],na.rm = T)/nb_aa_copie,
                              RSCU = sum(observation[codon]) / nb_observation
                            ))
          } else {
            
            aa_data = rbind(aa_data,
                            data.frame(
                              codon ,
                              nb_syn,
                              amino_acid,
                              rules = "watson-crick",
                              letter_aa = codon ,
                              anticodon = chartr("TUACG","AATGC",stri_reverse(codon)),
                              count = sum(observation[codon]) / sum(observation[codon_used],na.rm = T),
                              RTF= sum(tRNASE_copies_table[codon],na.rm = T)/nb_aa_copie,
                              RSCU = sum(observation[codon]) / nb_observation
                            ))
            
          }
        }
      }
      aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]
      aa_data = aa_data[aa_data$nb_syn >= 2 ,]
      
      aa_data_lm = aa_data[ aa_data$nb_syn != 2,] ## only TriQuaSex
      
      lm_model_obs = summary( lm( aa_data_lm$RTF ~ aa_data_lm$RSCU ) )
      
      p_value = lm_model_obs$coefficients[2,4]
      R2 = lm_model_obs$r.squared
      
      
      ######## EXPECTED CODON OPTIMAL         ########         ########         ########         ########         ########         ########         ######## 
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
      
      data.frame_graph = data.frame()
      
      length_limit="all"
      codon_usage_selected = codon_usage
      
      if (length_limit == "long"){  
        codon_usage_selected = codon_usage[codon_usage$length_cds >= median(codon_usage$length_cds),]
      } else if (length_limit == "short"){  
        codon_usage_selected = codon_usage[codon_usage$length_cds < median(codon_usage$length_cds),]
      }
      
      xaxis = codon_usage_selected$median_fpkm 
      proportion = 5/100
      quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
      q50 = median(xaxis,na.rm = T )
      intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
      
      FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
      # FPKM_bins = log10(FPKM_bins)
      
      GC3_obs = rowSums(codon_usage_selected[c("C3","G3")],na.rm = T)
      ATGC3_obs = rowSums(codon_usage_selected[c("A3","T3","C3","G3")],na.rm = T)
      
      GCi_obs = rowSums(codon_usage_selected[c("Ci","Gi")],na.rm = T)
      ATGCi_obs = rowSums(codon_usage_selected[c("Ai","Ti","Ci","Gi")],na.rm = T)
      
      data.frame_graph = rbind(data.frame_graph,data.frame(
        freq = tapply( GCi_obs / ATGCi_obs  , intervalle_FPKM , function(x) mean(x,na.rm=T)) ,
        fpkm = FPKM_bins,
        categorie = "GC intronic",
        group = "control",
        type_aa="GC3-GCi",
        set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
        nb_codon = NA,
        list=NA
      ))
      
      data.frame_graph = rbind(data.frame_graph,data.frame(
        freq = tapply( GC3_obs / ATGC3_obs  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
        fpkm = FPKM_bins,
        categorie = "GC3",
        group = "control",
        type_aa="GC3-GCi",
        set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
        nb_codon = NA,
        list= NA
      ))
      
      data.frame_graph = rbind(data.frame_graph,data.frame(
        freq = tapply( codon_usage_selected$length_cds  , intervalle_FPKM , function(x) median(x,na.rm=T)),
        fpkm = FPKM_bins,
        categorie = "CDS length",
        group = "control",
        type_aa="GC3-GCi",
        set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
        nb_codon = NA,
        list= NA
      ))
      # type_aa="TriQuaSex"
      # for ( type_aa in c("Wb_WC_notambiguous","WC_duet","WC_duet_ambiguous","WC_duet_notambiguous")){
      for ( type_aa in c("WB_vs_WC_notduet","IC","GU" , c("Wb_WC_notambiguous","WC_duet","WC_duet_ambiguous","WC_duet_notambiguous"))){
        
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
          
        } else if  (type_aa %in%c("IC","GU" )){
          dt_selected = code_table[ code_table$nb_syn > 2,]
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
          
        } else {
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
        
        
        data.frame_graph = rbind(data.frame_graph,data.frame(
          freq =  tapply( COA_obs / COA_neg_obs  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
          fpkm = FPKM_bins,
          categorie = "optimal codons",
          group = "not_control",
          type_aa=paste(type_aa," codons, Nb aa = ",nb_aa,sep=""),
          set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
          nb_codon ,
          list=paste(list_COA,collapse =";")
        ))
        
        
        data.frame_graph = rbind(data.frame_graph,data.frame(
          freq =   tapply( COA_obs_intronic / COA_neg_obs_intronic  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
          fpkm = FPKM_bins,
          categorie = "intronic triplets control",
          group = "control",
          type_aa=paste(type_aa," codons, Nb aa = ",nb_aa,sep=""),
          set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
          nb_codon ,
          list=paste(list_COA,collapse=";")
        ))
        
        
        
        value = tapply( ( COA_obs / COA_neg_obs )   , intervalle_FPKM , function(x) mean(x,na.rm=T)) -
          tapply(  ( COA_obs_intronic / COA_neg_obs_intronic ) , intervalle_FPKM , function(x) mean(x,na.rm=T))
        
        # value = tapply( ( COA_obs / COA_neg_obs )  - ( COA_obs_intronic / COA_neg_obs_intronic ) , intervalle_FPKM , function(x) mean(x,na.rm=T))
        value = ( value[length(value)] -  median(value[1:10])) / ( FPKM_bins[length(FPKM_bins)] -  median(FPKM_bins[1:10]))
        
        intervalle_FPKM
        xaxis < q50
        
        ecart_cds = (tapply( COA_obs / COA_neg_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)] - mean( COA_obs[xaxis < q50] / 
                                                                                                                                   COA_neg_obs[xaxis < q50] , na.rm=T))/
          mean( COA_obs[xaxis < q50] / 
                  COA_neg_obs[xaxis < q50] , na.rm=T)*100
        
        ecart_intron =(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)] - mean( COA_obs_intronic[xaxis < q50] / 
                                                                                                                                                       COA_neg_obs_intronic[xaxis < q50] , na.rm=T))/
          mean( COA_obs_intronic[xaxis < q50] / 
                  COA_neg_obs_intronic[xaxis < q50] , na.rm=T)*100
        
        
        somme = sum( COA_obs) / sum(COA_neg_obs ) - sum( COA_obs_intronic) / sum(COA_neg_obs_intronic )
        
        
        data.frame_graph = rbind(data.frame_graph,data.frame(
          freq =   tapply( ( COA_obs / COA_neg_obs )   , intervalle_FPKM , function(x) mean(x,na.rm=T)) -
            tapply(  ( COA_obs_intronic / COA_neg_obs_intronic ) , intervalle_FPKM , function(x) mean(x,na.rm=T)),
          fpkm = FPKM_bins,
          categorie = paste(type_aa,"slope :",round(value,4)),
          group = "Delta",
          type_aa = "Delta",
          set=paste(length_limit,", N = " ,vector_count[length_limit],sep=""),
          nb_codon ,
          list=paste(list_COA,collapse=";")
        ))
        
        gc3 = (GC3_obs / ATGC3_obs)
        gci = (GCi_obs / ATGCi_obs)
        model_lm = lm( gc3 ~  gci)
        
        spearman_method_gc3gci = cor.test( gc3, gci,method="spearman",exact=F)
        # print(spearman_method_gc3gci)
        
        data12 = rbind(data12,data.frame(
          species ,
          tRNA_GFF,
          # nb_busco = length(unique(busco_tab$busco_id)),
          rho_aa_fpkm = spearman_method_aa$estimate,
          pval_aa_fpkm = spearman_method_aa$p.value,
          nb_codon_not_decoded,
          slope = round(value,5),
          
          
          rho_gc3_gci = spearman_method_gc3gci$estimate,
          pval_gc3_gci = spearman_method_gc3gci$p.value,
          gc3_sum = sum(GC3_obs ,na.rm=T) / sum( ATGC3_obs,na.rm=T),
          gc3 = mean(GC3_obs / ATGC3_obs,na.rm=T),
          std_gc3 = stderror(GC3_obs / ATGC3_obs),
          var_gc3 = var(GC3_obs / ATGC3_obs,na.rm=T),
          
          gci_sum = sum(GCi_obs ,na.rm=T) / sum( ATGCi_obs,na.rm=T),
          gci = mean(GCi_obs / ATGCi_obs,na.rm=T),
          std_gci = stderror(GCi_obs / ATGCi_obs),
          var_gci = var(GCi_obs / ATGCi_obs,na.rm=T),
          
          gc3_top5 = tapply( GC3_obs / ATGC3_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
          var_gc3_top5 = tapply( GC3_obs / ATGC3_obs   , intervalle_FPKM , function(x) var(x,na.rm=T))[length(FPKM_bins)],
          gci_top5 = tapply( GCi_obs / ATGCi_obs   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],
          var_gci_top5 = tapply( GCi_obs / ATGCi_obs   , intervalle_FPKM , function(x) var(x,na.rm=T))[length(FPKM_bins)],
          
          lm_r2_gc3_gci = summary(model_lm)$r.squared,
          lm_pval_gc3_gci = summary(model_lm)$coefficients[2,4],
          
          average_RTF_abundant_tRNA = mean( code_table[list_COA,]$RTF[code_table[list_COA,]$RTF != 0] ) ,
          
          optifreq_top5 = round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5),
          opti_freq_low50 = round( mean( COA_obs[codon_usage_selected$median_fpkm <= median(codon_usage_selected$median_fpkm )] /
                                           COA_neg_obs[codon_usage_selected$median_fpkm <= median(codon_usage_selected$median_fpkm )] , na.rm=T) , 5),
          optifreq_top5_intron = round(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5),
          opti_freq_low50_intron = round( mean( COA_obs_intronic[codon_usage_selected$median_fpkm <= median(codon_usage_selected$median_fpkm )] /
                                                  COA_neg_obs_intronic[codon_usage_selected$median_fpkm <= median(codon_usage_selected$median_fpkm )] , na.rm=T),5),
          
          ecart_cds = round(ecart_cds,5),
          ecart_intron = round(ecart_intron,5),
          somme = round(somme,5),
          length_limit,
          type_aa,
          nb_genes = vector_count[length_limit],
          nb_aa,
          nb_scu=length(list_COA)
        ))
        
      }
    }
  }
}

write.table(data12 , paste("data/data12.tab",sep=""),sep="\t",quote=F,row.names=F)



