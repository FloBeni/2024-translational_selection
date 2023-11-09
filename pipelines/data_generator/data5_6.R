library(stringi)
library(ggplot2)
library(stringr)
library(ape)


path_data = "/home/fbenitiere/data/"

wobble_pairing = c("C"="IC","T"="GU","G"="UG","A"="IA")

wobble_associat_wc = list("T"="C",
                          "A"="A",
                          "G"="A",
                          "C"="T")

code = read.delim(paste(path_data,"Projet-SplicedVariants/Fichiers-data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")

amino_acid_list = unique(code[code$nb_syn > 1 & code$aa_name != "Ter",]$aa_name_scu)
amino_acid_list = unique(code[ code$aa_name != "Ter",]$aa_name)


# all_data6v1 = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75.tab",sep=""))
# all_data6v2 = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_prop_gap_pergene_25_50_75_rmfirst1000bp.tab",sep=""))
all_data6v1 = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_cons_pergene_25_50_75.tab",sep=""))
all_data6v2 = read.delim(paste(path_data,"Projet-NeGA/translational_selection/scu_on_constraint_site/compilation_cons_pergene_25_50_75_rmfirst1000bp.tab",sep=""))

species = "Homo_sapiens"


clade_dt = read.delim(paste( "data/clade_dt.tab",sep=""),header=T)
rownames(clade_dt) = clade_dt$species
clade_dt$clade_group = factor(clade_dt$clade_group, levels = c("Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Teleostei","Mammalia","Aves","Other Tetrapods"))


data12 = read.delim("data/data12.tab")
data12$clade_group = clade_dt[data12$species,]$clade_group
data12 = data12[ data12$nb_codon_not_decoded == 0 & data12$pval_aa_fpkm < 0.05 ,]
list_species = unique(data12$species)

data6 = data.frame()
data5 = data.frame()
for ( species in list_species){print(species)
  if (clade_dt[species,]$clade_group %in% c("Mammalia","Aves","Other Tetrapods")){
    all_data6_sub = all_data6v2[all_data6v2$species == species,] 
  } else {
    all_data6_sub = all_data6v1[all_data6v1$species == species,] 
  }
  
  codon_usage = read.delim( paste(path_data,"/Projet-SplicedVariants/Analyses/",species,"/codon_usage_gene_fpkm.tab",sep="") )
  codon_usage = codon_usage[codon_usage$protein_id %in% all_data6_sub$protein,]
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  
  codon_usage = codon_usage[codon_usage$intern_stop_codon==0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %%3 == 0,]
  
  if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){ 
    codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else {    print(species)}
  
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  
  
  all_data6_sub = all_data6_sub[all_data6_sub$protein %in% codon_usage$protein_id,]
  all_data6_sub$length =  (all_data6_sub$len_high_const_seq + all_data6_sub$len_mod_const_seq + all_data6_sub$len_slight_const_seq + all_data6_sub$len_unconst_seq)
  
  all_data6_sub = all_data6_sub[all_data6_sub$len_high_const_seq != 0 & all_data6_sub$len_unconst_seq != 0 , ]
  
  # print(ggplot(all_data6_sub) + geom_boxplot(aes(x="len_high_const_seq",y=len_high_const_seq/length))+
  #         geom_boxplot(aes(x="len_mod_const_seq",y=len_mod_const_seq/length))+
  #         geom_boxplot(aes(x="len_slight_const_seq",y=len_slight_const_seq/length))+
  #         geom_boxplot(aes(x="len_unconst_seq",y=len_unconst_seq/length))) +
  #   stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red",aes(x="len_high_const_seq",y=len_high_const_seq/length))+
  #   stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red",aes(x="len_mod_const_seq",y=len_mod_const_seq/length))+
  #   stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red",aes(x="len_slight_const_seq",y=len_slight_const_seq/length))+
  #   stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red",aes(x="len_unconst_seq",y=len_unconst_seq/length))
  
  code_table = read.delim(paste(path_data,"/Projet-NeGA/translational_selection/ExpOpti/",species,"_code_table.tab",sep=""))
  code_table$wb = wobble_pairing[substr(code_table$codon,3,3)]
  code_table$decoded_codon = sapply(code_table$codon,function(x) paste(substr(x,1,2),wobble_associat_wc[[substr(x,3,3)]],collapse="",sep=""))
  rownames(code_table) = code_table$codon
  code_table = code_table[code_table$aa_name %in% amino_acid_list,]
  nb_codon_not_decoded = sum(!code_table$decoded)
  
  type_aa = "Wb_WC_notambiguous"
  for ( type_aa in c("IC","WB_vs_WC_notduet","Wb_WC_notambiguous","WC_duet_ambiguous")){
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
      
    } else if  (type_aa %in% c("IC","GU" )){
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
      
    }  else {
      dt_selected = code_table[ code_table$aa_name == type_aa,]
      optimal_count = table( dt_selected[ dt_selected$Wobble_abond | dt_selected$WC_abond,]$aa_name )
      
      list_aa = names(optimal_count)[ optimal_count != table(code$aa_name)[names(optimal_count)]]
      nb_aa = length(list_aa)
      
      list_COA = dt_selected[(dt_selected$Wobble_abond | dt_selected$WC_abond) & dt_selected$aa_name %in% list_aa,]$codon
      list_COA_neg = code[code$aa_name %in% list_aa,]$codon
      nb_codon = length(list_COA)
      
    }
    
    
    if (length(list_COA) != 0){
      table_constrain = data.frame(busco_id = all_data6_sub$busco_id)
      for (constrain in c("_highconst","_modconst","_sligconst","_unconst")){
        table_constrain[,paste("COA",constrain,sep="")] = rowSums(all_data6_sub[paste(list_COA,constrain,sep = "")],na.rm = T)
        table_constrain[,paste("COA_neg",constrain,sep="")] = rowSums(all_data6_sub[paste(list_COA_neg,constrain,sep = "")],na.rm = T)
      }
      
      
      data5 = rbind(data5, data.frame(
        species,
        type_aa,
        freq = table_constrain$COA_highconst / table_constrain$COA_neg_highconst ,
        nb_site = sum(all_data6_sub$len_high_const_seq),
        nb_genes = nrow(codon_usage),
        categorie = "Highly constrained"
      ))
      
      data5 = rbind(data5, data.frame(
        species,
        type_aa,
        freq = table_constrain$COA_modconst / table_constrain$COA_neg_modconst ,
        nb_site = sum(all_data6_sub$len_mod_const_seq),
        nb_genes = nrow(codon_usage),
        categorie = "Moderately constrained"
      ))
      data5 = rbind(data5, data.frame(
        species,
        type_aa,
        freq = table_constrain$COA_sligconst / table_constrain$COA_neg_sligconst ,
        nb_site = sum(all_data6_sub$len_slight_const_seq),
        nb_genes = nrow(codon_usage),
        categorie = "Slighlty constrained"
      ))
      data5 = rbind(data5, data.frame(
        species,
        type_aa,
        freq = table_constrain$COA_unconst / table_constrain$COA_neg_unconst ,
        nb_site = sum(all_data6_sub$len_unconst_seq),
        nb_genes = nrow(codon_usage),
        categorie = "Unconstrained"
      ))
      
      tapply(data5$freq,data5$categorie,function(x) mean(x,na.rm=T))
      
      data6 = rbind(data6,data.frame(species,
                                     type_aa,
                                     nb_aa,
                                     nb_scu = length(list_COA),
                                     nb_genes = nrow(all_data6_sub),
                                     mean_highconstrain = mean(table_constrain$COA_highconst/table_constrain$COA_neg_highconst,na.rm = T),
                                     len_highconstrain = sum(table_constrain$COA_neg_highconst,na.rm = T),
                                     mean_modconstrain = mean(table_constrain$COA_modconst/table_constrain$COA_neg_modconst,na.rm = T),
                                     len_modconstrain = sum(table_constrain$COA_neg_modconst,na.rm = T),
                                     mean_sligconstrain = mean(table_constrain$COA_sligconst/table_constrain$COA_neg_sligconst,na.rm = T),
                                     len_sligconstrain = sum(table_constrain$COA_neg_sligconst,na.rm = T),
                                     mean_unconstrain = mean(table_constrain$COA_unconst/table_constrain$COA_neg_unconst,na.rm = T),
                                     len_unconstrain = sum(table_constrain$COA_neg_unconst,na.rm = T)
      ))
    }
  }
}
data5$categorie = factor(data5$categorie,levels = unique(data5$categorie))

data5 = data5[data5$species %in% c("Homo_sapiens","Caenorhabditis_elegans","Drosophila_melanogaster"),]

write.table(data5,"data/data5_cons.tab",quote=F,row.names = F,sep="\t")
write.table(data6,"data/data6_cons.tab",quote=F,row.names = F,sep="\t")
