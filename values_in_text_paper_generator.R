source("figure/figure_main_generator/library_path.R")

{
  data1 = read.delim("data/data1_supp.tab")
  rownames(data1) = data1$species
  data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
  w
  nrow(data1[ data1$nb_genes_filtered >= 5000,])
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
  print("Non-adaptive processes shape codon usage variations across species")
  print(paste("metazoans covering a wide range of clades (N=",nrow(data1),", ",sum(table(data1$clade_group)[c("Mammalia","Aves","Other Tetrapods","Teleostei")]),
              " Vertebrates and ",nrow(data1) - sum(table(data1$clade_group)[c("Mammalia","Aves","Other Tetrapods","Teleostei")])," Invertebrates)" , sep=""))
  
  print(paste("GCi extends from ",round(min(data1$gci),2)," to ",round(max(data1$gci),2),sep=""))
  print(paste("GC3 span from ",round(min(data1$gc3),2)," to ",round(max(data1$gc3),2),sep=""))
  
  
  print(paste("GC3 and GCi are highly correlated (rho=",round(data1["Homo_sapiens",]$rho_gc3_gci,2),")",sep=""))
  
  data2 = read.delim("data/data2_supp.tab")
  dt_graph = data2[data2$species == "Homo_sapiens" & !is.na(data2$GCi),]
  print(paste("(10th percentile ranging from ",round(quantile(dt_graph$GC3,0.1),2)," to ",round(quantile(dt_graph$GC3,0.9),2),")",sep=""))
  
  
  ##  Estimating the relative abundance of tRNA isodecoders
  print("Estimating the relative abundance of tRNA isodecoders")
  
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
  
  print(paste("the amino acid usage and the tRNA pool (rho=",round(data1["Caenorhabditis_elegans",]$rho_aa_fpkm,2),")",sep=""))
  
  print(paste("the studied species (N=",nrow(data1)," species)",sep=""))
  
  print(paste("in ",round(sum(data1$pval_aa_fpkm < 0.05)/nrow(data1)*100),"% of the analyzed species",sep=""))
  print(paste("good proxy of their abundance, for the rest of our analyses we are focusing on this set of species (N=",sum(data1$pval_aa_fpkm < 0.05),")",sep=""))
  
  
  ##  Identifying putative optimal codons for translation
  print("Identifying putative optimal codons for translation")
  
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
  
  
  data4 = read.delim("data/data4_supp.tab")
  data4 = data4[data4$species == "metazoa",]
  codon_count_absent = sapply(data4$codon,function(x){
    dt = data4[data4$codon == x & data4$Var1 %in% c("WBp + abond","WBp","not decoded"),]
    sum(dt$Freq)
  })
  tapply(data4$Freq,data4$codon,sum)[1]
  print(paste("these seven codons lack the corresponding tRNA in most species (average=",
              round(mean(codon_count_absent[c("AGT","TTT","AAT","GAT","CAT","TGT","TAT")]/tapply(data4$Freq,data4$codon,sum)[1] * 100)),"%; Fig 2B)",sep=""))
  
  print(paste("except for the Glycine where GGT tRNA is highly avoided (",round(mean(codon_count_absent[c("GGT")]/tapply(data4$Freq,data4$codon,sum)[1] * 100)),"%).",sep=""))
  
  
  data1 = data1[data1$pval_aa_fpkm < 0.05,]
  print(paste("On average, the optimal codons can be defined for ", round(mean(data1$nb_aa_WB_WC_notambiguous))," amino acids per species.",sep=""))
  print(paste("In our analysis, we excluded ",sum(data1$nb_codon_not_decoded != 0)," species in which certain",sep=""))
  
  
  ## Highly expressed genes are enriched in POC
  print("Highly expressed genes are enriched in POC")
  
  data1 = data1[ data1$nb_codon_not_decoded == 0 & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000,]
  print(paste("5,000 protein-coding genes expressed (N=",nrow(data1)," species)",sep=""))
  
  print(paste("The results depicted in Fig4C indicate that, for ",table(data1$expressed_overused_background_WB_WC_notambiguous > 0)["TRUE"],
              " species (",round(table(data1$expressed_overused_background_WB_WC_notambiguous > 0)["TRUE"]/nrow(data1)*100),"%)",sep=""))
  
  
  print(paste("Resulting in an elevation of ",round(data1[data1$species == "Caenorhabditis_elegans",]$expressed_overused_background_WB_WC_notambiguous),
              "% in Caenorhabditis elegans, ",round(mean(data1[data1$clade_group == "Mecopterida",]$expressed_overused_background_WB_WC_notambiguous)),
              "% in Mecopterida, and ",
              round(mean(data1[data1$clade_group %in% c("Invertebrates","Mammalia","Aves","Other Tetrapods","Teleostei"),]$expressed_overused_background_WB_WC_notambiguous)),
              "% in vertebrates, with Aves exhibiting a higher translational intensity.",sep=""))
  
  print(paste("Not accounting for POCMT variations does not changed the results as it remains positive for ",table(data1$expressed_overused_WB_WC_notambiguous > 0)["TRUE"]," species ",sep=""))
  
  
  ## The most constraint sites tend to be enriched in POC
  print("The most constraint sites tend to be enriched in POC")
  
  data6 = read.delim("data/data6_supp.tab")
  dt_analysis = data6[data6$species == "Caenorhabditis_elegans",]
  
  print(paste("the most constrained sites displaying a higher proportion of POC (on average ",
              round(mean(dt_analysis[dt_analysis$categorie == "Highly constrained",]$freq,na.rm=T)*100,1),"% vs ",
              round(mean(dt_analysis[dt_analysis$categorie == "Unconstrained",]$freq,na.rm=T)*100,1),"%;",sep=""))
  
  print(paste("In Mecopterida we observe a difference reaching ",round(max(data1[data1$clade_group == "Mecopterida",]$constraint_overused_WB_WC_notambiguous),1),
              "% (with an average of ",round(mean(data1[data1$clade_group == "Mecopterida",]$constraint_overused_WB_WC_notambiguous),1),
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
