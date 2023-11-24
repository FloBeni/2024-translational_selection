source("figure/figure_main_generator/library_path.R")

{
  
  ##  Introduction
  data1 = read.delim("data/data1.tab")
  data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
  
  print(paste("selection intensity across " , nrow(data1) , " metazoan species" , sep = ""))
  
  
  ##  Estimating the relative abundance of tRNA isodecoders
  print(paste("metazoans covering a wide range of clades (N=",nrow(data1),", ",sum(table(data1$clade_group)[c("Mammalia","Aves","Other Tetrapods","Teleostei")]),
              " Vertebrates and ",nrow(data1) - sum(table(data1$clade_group)[c("Mammalia","Aves","Other Tetrapods","Teleostei")])," Invertebrates)" , sep=""))
  
  
  print(paste("but in some cases (",table(dt_graph$tRNA_GFF)["not from GFF"]," species over ",nrow(dt_graph),
              ") the tRNA are not annotated. We annotated the remaining ",table(dt_graph$tRNA_GFF)["not from GFF"],sep=""))
  
  
  dt_graph = dt_graph[ dt_graph$nb_gene > 5000 , ]
  print(paste("the studied species with at least 5,000 protein-coding genes expressed (N=",nrow(dt_graph),")",sep=""))
  
  
  print(paste("In ",round(sum(dt_graph$pval_aa_fpkm < 0.05)/nrow(dt_graph)*100),"% of the analyzed species",sep=""))
  print(paste("good proxy of their abundance, for the rest of our analyses we are focusing on this set of species (N=",sum(dt_graph$pval_aa_fpkm < 0.05),")",sep=""))
  
  ## Identifying codons optimized for translation (optimal codons)
  
  species = "Homo_sapiens"
  path = "/home/fbenitiere/data/"
  
  code = read.delim(paste(path,"Projet-SplicedVariants/Fichiers-data/standard_genetic_code.tab",sep=""))
  rownames(code) = code$codon
  code$nb_syn = table(code$aa_name)[code$aa_name]
  code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
  
  codon_usage = read.delim( paste(path,"/Projet-SplicedVariants/Analyses/",species,"/codon_usage_gene_fpkm.tab",sep="") )
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  codon_usage = codon_usage[codon_usage$length_cds == codon_usage$length,]
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage$stop_codon = rowSums(codon_usage[,stop_codon])
  codon_usage = codon_usage[codon_usage$stop_codon == 1,]
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 & !is.na(codon_usage$median_fpkm) ,]
  
  
  observation = colSums( codon_usage[3:70] , na.rm = T )
  AAC = observation["AAC"]
  AAT = observation["AAT"]
  print(paste("AAT accounting for ",round(AAT/(AAC+AAT)*100),"% of the occurrences, respectively.",sep=""))
  
  
  
  data2 = read.delim("data/data2.tab")
  data2 = data2[data2$species == "metazoa",]
  codon_count_absent = sapply(data2$codon,function(x){
    dt = data2[data2$codon == x & data2$Var1 %in% c("WB + abond","WB","not decoded"),]
    sum(dt$Freq)
  })
  tapply(data2$Freq,data2$codon,sum)[1]
  print(paste("For these 7 codons the corresponding tRNA is absent in most species (average = ",round(mean(codon_count_absent[c("TAT","TGT","CAT","AAT","GAT","TTT","AGT")]/tapply(data2$Freq,data2$codon,sum)[1] * 100)),"%; Fig 2B)",sep=""))
  
  
  print(paste("except for the Glycine where GGT tRNA is highly avoided (",round(mean(codon_count_absent[c("GGT")]/tapply(data2$Freq,data2$codon,sum)[1] * 100)),"%).",sep=""))
  
  
  Freq_opti = read.delim("data/Freq_opti.tab")
  Freq_opti = Freq_opti[Freq_opti$type_aa == "Wb_WC_notambiguous",]

  
  arbrePhylo = read.tree(paste("data/phylogenetic_tree_root.nwk",sep=""))
  list_species = arbrePhylo$tip.label
  Freq_opti = Freq_opti[Freq_opti$species %in% dt_graph$species[dt_graph$pval_aa_fpkm < 0.05],]  
  table(Freq_opti$tRNA_GFF)
 
  print(paste("On average, the optimal codons can be defined for ", round(mean(Freq_opti$nb_aa))," amino acids per species.",sep=""))
  print(paste("In our analysis, we excluded ",sum(Freq_opti$nb_codon_not_decoded != 0)," species in which certain",sep=""))
  
  
  
  
  
  ## Highly expressed genes are enriched in optimal codons
  
  Freq_opti = Freq_opti[Freq_opti$nb_codon_not_decoded == 0 ,]
  
  Freq_opti$ecart = (Freq_opti$optifreq_top5-Freq_opti$opti_freq_low50) - (Freq_opti$optifreq_top5_intron-Freq_opti$opti_freq_low50_intron)
  
  print(paste("The results depicted in Fig3E indicate that, for ",table(Freq_opti$ecart > 0)["TRUE"]," species (",round(table(Freq_opti$ecart > 0)["TRUE"]/nrow(Freq_opti)*100),"%)",sep=""))
  
  Freq_opti$clade_group = clade_dt[Freq_opti$species,]$clade_group
  
  
  print(paste("Resulting in an elevation of ",round(Freq_opti[Freq_opti$species == "Caenorhabditis_elegans",]$ecart*100),
              "% in Caenorhabditis elegans, ",round(mean(Freq_opti[Freq_opti$clade_group == "Lepido Diptera",]$ecart*100)),
              "% in Lepido Diptera, and ",round(mean(Freq_opti[Freq_opti$clade_group %in% c("Invertebrates","Mammalia","Aves","Other Tetrapods","Teleostei"),]$ecart*100)),
              "% in vertebrates, with Aves exhibiting a higher translational intensity.",sep=""))
  
  Freq_opti$ecart = (Freq_opti$optifreq_top5-Freq_opti$opti_freq_low50)
  
  print(paste("Accounting for OCMT variations do not change the results, positive for ",table(Freq_opti$ecart > 0)["TRUE"]," species",sep=""))
  
  
  ## The most constraint sites tend to be enriched in optimal codons
  
  data6 = read.delim("data/data6_gap33_pergene.tab")
  dt_analysis = data6[data6$type_aa == "Wb_WC_notambiguous" & data6$species == "Caenorhabditis_elegans",]
  print(paste("the most constrained sites displaying a higher proportion of OC (on average ",round(dt_analysis$mean_modconstrain*100,1),"% vs ",round(dt_analysis$mean_unconstrain*100,1),"%;",sep=""))
  
  data6 = data6[data6$species %in% Freq_opti$species,]
  
  dt_analysis = data6[data6$type_aa == "Wb_WC_notambiguous" ,]
  dt_analysis$ecart = dt_analysis$mean_modconstrain - dt_analysis$mean_unconstrain
  
  dt_analysis$clade_group = clade_dt[dt_analysis$species,]$clade_group
  print(paste("For most of the species (",sum(dt_analysis$ecart>0)/nrow(dt_analysis)*100,"%) the shift is positive suggesting that constraint sites",sep=""))
  
  
  print(paste("We observe that in Diptera there is a differences reaching ",round(max(dt_analysis[dt_analysis$clade_group == "Lepido Diptera",]$ecart*100)),
              "% (average = ",round(mean(dt_analysis[dt_analysis$clade_group == "Lepido Diptera",]$ecart*100)),"%.",sep=""))
  
  
}










