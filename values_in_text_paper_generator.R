source("figure/figure_main_generator/library_path.R")

{
  data1 = read.delim("data/data1_supp.tab")
  rownames(data1) = data1$species
  data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
  
  nrow(data1[  data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000,])
  nrow(data1[ data1$nb_codon_not_decoded == 0 & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000,])
  
  ##  Abstract
  print("Abstract")
  print(paste("encompassing " , nrow(data1) , " species" , sep = ""))
  print(paste("utilizing up to ",round(max(data1$expressed_overused_background_global)),"% more optimal codons" , sep = ""))
  
  
  ##  Introduction
  print("Introduction")
  print(paste("selection intensity across " , nrow(data1) , " metazoan species" , sep = ""))
  
  
  ##  Non-adaptive processes shape codon usage variations across species
  print("Non-adaptive processes are the primary drivers of codon usage variations among metazoans")
  
  print(paste("position between ",nrow(data1[ data1$nb_genes_filtered >= 5000,])," metazoan species",sep=""))
  print(paste("after removal of species lacking gene expression data (N=",nrow(data1[ data1$nb_genes_filtered < 5000,]),")",sep=""))
  
  data1 = data1[ data1$nb_genes_filtered >= 5000,]
  
  print(paste("clades (",sum(table(data1$clade_group)[c("Mammalia","Aves","Other Tetrapods","Teleostei")]),
              " vertebrates, ",sum(table(data1$clade_group)[c("Diptera","Lepidoptera","Coleoptera","Hymenoptera","Other Insects")]),
              " insects and ",sum(table(data1$clade_group)[c("Other Metazoans","Nematoda")]),
              " other metazoans).",sep=""))
  
  
  
  
  print(paste("GCi extends from ",round(min(data1$gci),2)," to ",round(max(data1$gci),2),sep=""))
  print(paste("GC3 span from ",round(min(data1$gc3),2)," to ",round(max(data1$gc3),2),sep=""))
  
  
  
  print(paste("Diptera clade (N=",sum(table(data1$clade_group)[c("Diptera")]),")",sep=""))
  print(paste("the wider range of GC3 variations (from ",round(min(data1[data1$clade_group == "Diptera",]$gc3),2),"% to ",round(max(data1[data1$clade_group == "Diptera",]$gc3),2),"%).",sep=""))
  
  
  print(paste("GC3 and GCi are highly correlated (rho=",round(data1["Homo_sapiens",]$rho_gc3_gci,2),")",sep=""))
  print(paste("pronounced in Caenorhabditis elegans (rho=",round(data1["Caenorhabditis_elegans",]$rho_gc3_gci,2),",",sep=""))
  
  
  ##  tRNA abundance matches transcriptome requirements
  print("tRNA abundance matches transcriptome requirements")
  
  print(paste("but in some cases (",sum(!data1$tRNA_GFF)," species over ",nrow(data1),
              ") the tRNA are not annotated. We annotated the remaining ",sum(!data1$tRNA_GFF),sep=""))
  
  data7 = read.delim("data/data7_supp.tab")
  dt_graph = data7[data7$species == "Drosophila_melanogaster",]
  spearman_method_aa = cor.test( dt_graph$prop_abundance_average, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
  print(paste(" (rho = ",round(spearman_method_aa$estimate,2),")",sep=""))
  dt_graph = data7[data7$species == "Homo_sapiens",]
  spearman_method_aa_human = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
  dt_graph = data7[data7$species == "Drosophila_melanogaster",]
  spearman_method_aa_droso = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
  print(paste(" (rho = ",round(spearman_method_aa$estimate,2),") as well with the gene copy number, (rho = ",round(spearman_method_aa_droso$estimate,2)," and ",round(spearman_method_aa_human$estimate,2)," respectively)",sep=""))
  
  print(paste(" This analysis was conducted across all studied species (N=",nrow(data1)," species )",sep=""))
  
  print(paste("in ",round(sum(data1$pval_aa_fpkm < 0.05)/nrow(data1)*100),"% of the analyzed species",sep=""))
  
  print(paste("with amino acid usage (N=",sum(data1$pval_aa_fpkm < 0.05)," species)",sep=""))
  
  data1 = data1[data1$pval_aa_fpkm < 0.05,]
  
  
  ##  tRNA abundance defines putative-optimal codons
  print("tRNA abundance defines putative-optimal codons")
  
  species = "Homo_sapiens"
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
  rownames(code) = code$codon
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]
  codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),]
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
  
  observation = colSums( codon_usage[3:70] , na.rm = T )
  AAC = observation["AAC"]
  AAT = observation["AAT"]
  print(paste("AAT accounting for ",round(AAT/(AAC+AAT)*100),"% of the occurrences, respectively.",sep=""))
  
  print(paste("we excluded ",sum(data1$nb_codon_not_decoded != 0)," species in which certain",sep=""))
  data1 =   data1[ data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 & data1$nb_codon_not_decoded == 0,]
  
  data12 = read.delim("data/data12_supp.tab")
  data12 = data12[data12$species %in% data1$species,]
  
  list_species = c()
  for(species in unique(data12$species)){
    aa_present = table(data12[data12$species == species & data12$nb_syn != 1 & data12$abundance != 0 ,]$amino_acid)
    aa_present = names(aa_present[aa_present >= 2 ])
    if (any(!(aa_present %in% unique(data12[data12$species == species & data12$amino_acid %in% aa_present & data12$POC1, ]$amino_acid)))){
      print(species)
      list_species=append(list_species,species)
      print(aa_present[!(aa_present %in% unique(data12[data12$species == species & data12$amino_acid %in% aa_present & data12$POC1, ]$amino_acid))])
    }
  }
  
  print(paste(" gene copy number and decodes all synonymous codons (found in N=",length(list_species)," species)",sep=""))
  
  print(paste(" Phe, Asn, Asp, His, Cys and Tyr ( ",100*round(1-1/sum(table(data12[data12$POC2 & data12$species%in%data1$species,]$amino_acid)),3) ,")%",sep=""))
  
  
  print(paste("In the case of human we are able to identify POC1 for ",length(unique(data12[data12$species == "Homo_sapiens" & data12$POC1,]$amino_acid)),
              " amino acid and POC2 for ",length(unique(data12[data12$species == "Homo_sapiens" & data12$POC2,]$amino_acid))," codons.",sep=""))
  
  
  
  print(paste("we can define POC1 and POC2 for ",length(unique(data12[data12$species == "Caenorhabditis_elegans" & data12$POC1,]$amino_acid)),
              " and ",length(unique(data12[data12$species == "Caenorhabditis_elegans" & data12$POC2,]$amino_acid))," amino acids",sep=""))
  
  
  print(paste("POC1 are defined for ", round(mean(data1$nb_aa_POC1),1),
              " amino acids \textit{per} species (ranging from ",min(data1$nb_aa_POC1)," to ",max(data1$nb_aa_POC1),") and POC2 for ", round(mean(data1$nb_aa_POC2),1),
              " amino acids (ranging from ",min(data1$nb_aa_POC2)," to ",max(data1$nb_aa_POC2),")"
              ,sep="")) # Tyto alba is the only one having 7 because Ile
  
  ## Highly expressed genes are enriched in POC1 and POC2
  print("Highly expressed genes are enriched in POC1 and POC2")
  
  data5 = read.delim("data/data5_supp.tab")
  data5$gene_set = str_replace_all(data5$gene_set,"all,","genes,")
  data5[data5$categorie == "POC-matching triplets (POCMT)",]$categorie = "control"
  data5[data5$categorie == "Putative optimal codons (POC)",]$categorie = ""
  
  dt_graph = data5[ data5$species == "Caenorhabditis_elegans",]
  
  print(paste("we observed a rise in POC1 from ",round(min(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) &
                                                                      !grepl("control",dt_graph$categorie) &
                                                                      grepl("POC1",dt_graph$set),]$freq),2)*100,"% to ",
              round(dt_graph[!grepl("control",dt_graph$categorie) & dt_graph$fpkm == max(dt_graph$fpkm) & grepl("POC1",dt_graph$set) , ]$freq,2)*100,"% for highly expressed genes"))
  
  print(paste("This analysis was performed for each of the studied species (N=",nrow(data1)," species).",sep=""))
  
  print(paste("For ",table(data1$expressed_overused_background_POCs > 0)["TRUE"]," species (",round(table(data1$expressed_overused_background_POCs > 0)["TRUE"]/nrow(data1)*100),"%), the prevalence of POCs",sep=""))
  
  print(paste("The strongest variation is observed in Caenorhabditis elegans (+",
              round(data1[data1$species == "Caenorhabditis_elegans",]$expressed_overused_background_POCs),
              "%).",sep=""))
  
  
  print(paste("clades from +",round(mean(data1[data1$clade_group == "Diptera",]$expressed_overused_background_POCs)),
              "% on average in Diptera to and +",
              round(mean(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$expressed_overused_background_POCs)),
              "% in vertebrates",sep=""))
  
  print(paste("Not accounting for POC-control variations does not changed the results as it remains positive for ",table(data1$expressed_overused_POCs > 0)["TRUE"]," species ",sep=""))
  
  
  # print(paste(,sep=""))
  
  
  ## Most constrained sites tend to be enriched in POCs
  print("Most constrained sites tend to be enriched in POCs")
  
  
  # print(paste("for most tetrapods (N=",nrow(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods"),])," species)",sep=""))
  
  data6 = read.delim("data/data6_supp.tab")
  dt_analysis = data6[data6$species == "Caenorhabditis_elegans",]
  
  print(paste("the most constrained sites displaying a higher proportion of POCs (on average ",
              round(mean(dt_analysis[dt_analysis$type_aa == "POCs" & dt_analysis$categorie == "Highly constrained",]$freq,na.rm=T)*100,1),"% vs ",
              round(mean(dt_analysis[dt_analysis$type_aa == "POCs" & dt_analysis$categorie == "Unconstrained",]$freq,na.rm=T)*100,1),"%;",sep=""))
  
  print(paste("Our findings indicate that a majority of the species (",round(sum(data1$constraint_overused_POCs>0)/nrow(data1)*100,),"%) exhibit a positive shift in POC",sep=""))
  
  
  print(paste("In Diptera we observed a difference reaching ",round(max(data1[data1$clade_group == "Diptera",]$constraint_overused_POCs),1),
              "% (with an average of ",round(mean(data1[data1$clade_group == "Diptera",]$constraint_overused_POCs),1),
              "%).",sep=""))
  
  
  mean(dt_analysis[dt_analysis$categorie == "Highly constrained",]$freq,na.rm=T)*100 - mean(dt_analysis[dt_analysis$categorie == "Unconstrained",]$freq,na.rm=T)*100 
  
  data1[data1$species == "Caenorhabditis_elegans",]$constraint_overused_WB_WC_notambiguous
  
  
  ## Expressed genes codon usage modeled by tRNA pool
  print("Expressed genes codon usage modeled by tRNA pool")
  print(paste("Mecopterida exhibit on average ",round(mean(data1[data1$clade_group == "Mecopterida",]$expressed_overused_background_WB_WC_notambiguous)),
              "% difference in POC",sep=""))
  
  data4 = read.delim("data/data4_supp.tab")
  dt_graph = data4[ data4$method_to_calculate == "per_gene" & data4$group == "Wb_WC_notambiguous" ,]
  print(paste("for high expressed genes, from ",round(min(dt_graph$density_nonoptimal_to_optimal_codon*100),1),"% to ",
              round(max(dt_graph$density_nonoptimal_to_optimal_codon*100),1),"%. ",sep=""))
  
  dt_graph = data4[ data4$method_to_calculate == "per_gene" & data4$group == "duet_ambiguous" ,]
  print(paste("to wobble increases from ",round(min(dt_graph$density_nonoptimal_to_optimal_codon*100),1),"% to ",
              round(max(dt_graph$density_nonoptimal_to_optimal_codon*100),1),"%, suggesting that WBp to WCp",sep=""))
  
  dt_graph = data4[ data4$method_to_calculate == "per_gene" & data4$group == "ic_abondant" ,]
  print(paste("were preferentially selected in highly expressed genes (from ",round(min(dt_graph$density_optimal_to_nonoptimal_codon*100),1),"% to ",
              round(max(dt_graph$density_optimal_to_nonoptimal_codon*100),1),"%.",sep=""))
  
  
  ## tRNA pool partially modulated by neutral substitutions pattern
  print("tRNA pool partially modulated by neutral substitutions pattern")
  
  print(paste("difference in non-adaptive exposition (ranging from ",round(min(data1[data1$clade_group == "Mecopterida",]$gci),2),
              " to ",round(max(data1[data1$clade_group == "Mecopterida",]$gci),2),", ",sep=""))
  
  
  ## Discussion
  print("Discussion")
  
  data1 = read.delim("data/data1_supp.tab")
  print(paste("encompassing ",nrow(data1)," species.",sep=""))
  
}
