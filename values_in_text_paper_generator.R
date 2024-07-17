{
  source("figure/figure_main_generator/library_path.R")
  
  
  data1 = read.delim("data/data1_supp.tab",comment.char = "#")
  rownames(data1) = data1$species
  data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
  
  
  code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
  rownames(code) = code$codon
  code = code[code$aa_name != "Ter",]
  
  ##  Non-adaptive processes are the primary drivers of codon usage variations among metazoans
  print("Non-adaptive processes are the primary drivers of codon usage variations among metazoans")
  
  print(paste("initially selected ",nrow(data1)," metazoan species available",sep=""))
  
  
  print(paste("but we excluded ",sum(data1$nb_genes_filtered<5000)," species for which there were not enough transcriptomic data",sep=""))
  
  
  print(paste("base composition in the ",sum(data1$nb_genes_filtered>=5000)," remaining species",sep=""))
  
  data1 = data1[data1$nb_genes_filtered>=5000,]
  
  print(paste(sum(data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"))," vertebrates, ",
              sum(data1$clade_group %in% c("Diptera","Lepidoptera","Coleoptera","Hymenoptera","Other Insects")),
              " insects and ",sum(data1$clade_group %in%c("Other Metazoans","Nematoda"))," other metazoan species",sep=""))
  
  
  print(paste("Spearman's correlation coefficient, rho=",round(data1["Homo_sapiens",]$rho_gc3_gci,2),", p",sep=""))
  
  print(paste("less pronounced in Caenorhabditis elegans (rho=",round(data1["Caenorhabditis_elegans",]$rho_gc3_gci,2),",",sep=""))
  
  
  print(paste("in tetrapods (",sum(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods"),]$rho_gc3_gci > 0.7),"/",
              nrow(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods"),])," species with rho>0.7",
              sep=""))
  
  
  print(paste("in hymenopterans (",sum(data1[data1$clade_group %in% c("Hymenoptera"),]$rho_gc3_gci > 0.7),"/",
              nrow(data1[data1$clade_group %in% c("Hymenoptera"),])," species with rho>0.7",
              sep=""))
  
  
  print(paste("But overall, ",sum(data1$pval_gc3_gci < 0.05),"/",
              nrow(data1)," species (",sum(data1$pval_gc3_gci < 0.05)/nrow(data1)*100,
              "%", sep=""))
  
  ##  tRNA abundance matches proteome requirements
  print("tRNA abundance matches proteome requirements")
  
  data7 = read.delim("data/data7_supp.tab",comment.char = "#")
  data7 = data7[,code$anticodon]
  data7$clade_group = GTDrift_list_species[rownames(data7),]$clade_group
  
  print(paste("ranging from an average of ",
              round(mean(rowSums(data7[data7$clade_group == 'Hymenoptera',code$anticodon])),0),
              " tRNA gene copies per genome in hymenopterans to ",
              format(round(mean(rowSums(data7[data7$clade_group == 'Teleostei',code$anticodon])),0),big.mark=",",scientific=F)
              ," copies in teleost fish",sep=""))
  
  data7 = data7[,-which(colnames(data7)=="clade_group")]
  print(paste("Blattella germanica} contains ",data7[ "Blattella_germanica","AGA"]," copies of the AGA Ser-tRNA, vs. ",min(data7[ "Blattella_germanica",][data7[ "Blattella_germanica",] != 0 & data7[ "Blattella_germanica",] != max(data7[ "Blattella_germanica",])]),
              " to ",
              max(data7[ "Blattella_germanica",][data7[ "Blattella_germanica",] != 0 & data7[ "Blattella_germanica",] != max(data7[ "Blattella_germanica",])]),
              " copies for the other tRNAs",sep=""))
  
  
  data10 = read.delim("data/data10_supp.tab",comment.char = "#")
  dt_graph = data10[data10$species == "Drosophila_melanogaster",]
  spearman_method_aa = cor.test( dt_graph$prop_abundance_average, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
  print(paste("direct measures of tRNA abundance (rho = ",
              round(spearman_method_aa$estimate,2),sep=""))
  dt_graph = data10[data10$species == "Homo_sapiens",]
  spearman_method_aa_human = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
  dt_graph = data10[data10$species == "Drosophila_melanogaster",]
  spearman_method_aa_droso = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
  print(paste(" with their tRNA gene copy numbers (rho=",round(spearman_method_aa_droso$estimate,2)," and ",round(spearman_method_aa_human$estimate,2)," respectively",sep=""))
  
  
  print(paste(" same analysis conducted across ",nrow(data1)," animal species ",sep=""))
  
  print(paste("in ",round(sum(data1$pval_aa_fpkm<=0.05)/nrow(data1),2)*100,"% of the species ",sep=""))
  
  print(paste("correlates with amino acid usage (N=",sum(data1$pval_aa_fpkm<=0.05)," species)",sep=""))
  
  data1 = data1[data1$pval_aa_fpkm<=0.05,]
  print(paste("also excluded ",sum(data1$nb_codon_not_decoded!=0)," species ",sep=""))
  
  data1 = data1[data1$nb_codon_not_decoded == 0,]
  
  ##  Definition of putative-optimal codons based on tRNA abundance and wobble-pairing rules
  print("Definition of putative-optimal codons based on tRNA abundance and wobble-pairing rules")
  
  print(paste(nrow(code)," codons to their cognate tRNA",sep=""))
  
  data7 = data7[rownames(data7)%in%data1$species,]
  print(paste("ranges from ",min(rowSums(data7[,code$anticodon] != 0)),
              " to ",max(rowSums(data7[,code$anticodon] != 0)),
              " per species (average=",round(mean(rowSums(data7[,code$anticodon] != 0)),0),
              ").",sep=""))
  
  print(paste("This implies that ",61-max(rowSums(data7[,code$anticodon] != 0)),
              " to ",61-min(rowSums(data7[,code$anticodon] != 0)),sep=""))
  
  
  species = "Homo_sapiens"
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
  rownames(code) = code$codon
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.txt.gz",sep=""),comment.char = "#")
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
  print(paste("AAT accounts for ",round(AAT/(AAC+AAT)*100),"% of asparagine codons, highlighting the significance of wobble pairing",sep=""))
  
  print(paste("There are ",length(unique(code[code$nb_syn>1 & code$aa_name != "Ter",]$aa_name))," amino acids that",sep=""))
  
  
  
  data4 = read.delim("data/data4_supp.tab",comment.char = "#")
  data4 = data4[data4$species %in% data1$species,]
  
  list_species = c()
  list_exaequo = c()
  list_aa_present = c()
  for(species in unique(data4$species)){
    aa_present = table(data4[data4$species == species & data4$nb_syn != 1 & data4$tRNA_gene_copy != 0 ,]$amino_acid)
    aa_present = names(aa_present[aa_present >= 2 ])
    if (any(!(aa_present %in% unique(data4[data4$species == species & data4$amino_acid %in% aa_present & data4$POC1, ]$amino_acid)))){
      list_species=append(list_species,species)
      list_aa_present=append(list_aa_present,aa_present)
      list_exaequo=append(list_exaequo,aa_present[!(aa_present %in% unique(data4[data4$species == species & data4$amino_acid %in% aa_present & data4$POC1, ]$amino_acid))])
    }
  }
  
  print(paste("(",length(list_exaequo)," ex æquo among ",length(list_aa_present)," cases in total, in ",length(list_species)," species)",sep=""))
  
  print(paste("For the human genome, POC1 have been defined for ",length(unique(data4[data4$species == "Homo_sapiens" & data4$POC1,]$amino_acid)),
              " amino acids and POC2 for ",length(unique(data4[data4$species == "Homo_sapiens" & data4$POC2,]$amino_acid)),
              " amino acids",sep=""))
  
  
  print(paste("In contrast, for Caenorhabditis elegans, POC1 and POC2 are defined for ",length(unique(data4[data4$species == "Caenorhabditis_elegans" & data4$POC1,]$amino_acid)),
              " and ",length(unique(data4[data4$species == "Caenorhabditis_elegans" & data4$POC2,]$amino_acid)),
              " amino acids",sep=""))
  
  
  print(paste("On average among the ",nrow(data1),
              " species, POC1 are defined for ", round(mean(data1$nb_aa_POC1),1),
              " amino acids per species (ranging from ",min(data1$nb_aa_POC1)," to ",max(data1$nb_aa_POC1),
              ") and POC2 for ", round(mean(data1$nb_aa_POC2),1),
              " amino acids (ranging from ",min(data1$nb_aa_POC2)," to ",max(data1[data1$species!="Tyto alba",]$nb_aa_POC2)-1,
              ", except Tyto alba with ",data1["Tyto_alba",]$nb_aa_POC2," including Ile)",sep=""))
  
  ##  Highly expressed genes are enriched in optimal codons
  print("Highly expressed genes are enriched in optimal codons")
  
  data5 = read.delim("data/data5_supp.tab",comment.char = "#")
  
  dt_graph = data5[ data5$species == "Caenorhabditis_elegans",]
  
  print(paste("both for POC1 (from ",
              round(min(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) &
                                   grepl("codon",dt_graph$categorie) &
                                   grepl("POC1",dt_graph$set),]$freq),2)*100,"% to ",
              round(dt_graph[grepl("codon",dt_graph$categorie) & 
                               dt_graph$fpkm == max(dt_graph$fpkm) & 
                               grepl("POC1",dt_graph$set) , ]$freq,2)*100,"%) and for POC2 (from ",
              round(min(dt_graph[dt_graph$fpkm <= median(dt_graph$fpkm) &
                                   grepl("codon",dt_graph$categorie) &
                                   grepl("POC2",dt_graph$set),]$freq),2)*100,"% to ",
              round(dt_graph[grepl("codon",dt_graph$categorie) & 
                               dt_graph$fpkm == max(dt_graph$fpkm) & 
                               grepl("POC2",dt_graph$set) , ]$freq,2)*100,"% )",sep=""))
  
  
  
  
  species = "Caenorhabditis_elegans"
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
  rownames(code) = code$codon
  stop_codon = rownames(code[code$aa_name == "Ter",])
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.txt.gz",sep=""),comment.char = "#")
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]
  codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),]
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
  
  
  print(paste("which represent ",round(sum(codon_usage$median_fpkm <= 50)/nrow(codon_usage),2)*100,"% of genes in \textit{C. elegans}",sep=""))
  
  
  print(paste("the studied species (N=",nrow(data1)," species)",sep=""))
  
  summ_lm = summary(lm(data1$expressed_overused_background_POC1~data1$expressed_overused_background_POC2))
  print(paste("are strongly correlated (R2=",round(summ_lm$r.squared,2)*100,"% , p-value<",
              format(summ_lm$coefficients[2,4],scientific = T),sep=""))
  
  print(paste(" For ",sum(data1$expressed_overused_background_POC1 >=0,na.rm = T )," species (",
              round(sum(data1$expressed_overused_background_POC1 >=0,na.rm = T )/nrow(data1),2)*100,"%),",sep=""))
  
  
  print(paste(" For ",sum(data1$expressed_overused_background_POC2 >=0,na.rm = T )," species (",
              round(sum(data1$expressed_overused_background_POC2 >=0,na.rm = T )/nrow(data1),2)*100,"%;",sep=""))
  
  print(paste("are around +",round(mean(data1[data1$clade_group == "Diptera",]$expressed_overused_background_POCs)),
              "% in Diptera compared to +",
              round(mean(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$expressed_overused_background_POCs)),
              "% in vertebrates",sep=""))
  
  
  
  
  ##  Weak relationship between the strength of translational selection and the effective population size
  print("Weak relationship between the strength of translational selection and the effective population size")
  
  # Load data from Lynch et al. 2023
  # from https://www.embopress.org/doi/suppl/10.15252/embr.202357561/suppl_file/embr202357561-sup-0002-datasetev1.xlsx to which was added C.nigoni Ne
  lynch_dt = read.table("data/Lynch2023_embr202357561-sup-0002-metazoa.csv",header=T,sep="\t",dec=",")
  lynch_dt$species = str_replace_all(lynch_dt$Species," ","_")
  lynch_dt$species = sapply(lynch_dt$species ,function(x) paste(str_split_1(x,"_")[1],str_split_1(x,"_")[2],sep="_"))
  lynch_dt$genus = sapply(lynch_dt$species ,function(x) str_split_1(x,"_")[1])
  rownames(lynch_dt) = lynch_dt$species
  Ne_genus = tapply(lynch_dt$Ne,lynch_dt$genus,mean)
  
  data1$Ne = lynch_dt[data1$species,]$Ne
  data1$Ne_estimate = "from genus"
  data1[!is.na(data1$Ne),]$Ne_estimate = "from species"
  data1[is.na(data1$Ne),]$Ne = Ne_genus[sapply(data1[is.na(data1$Ne),]$species ,function(x) str_split_1(x,"_")[1])]
  data1[is.na(data1$Ne),]$Ne_estimate = ""
  
  
  print(paste("This list included ",sum(data1$Ne_estimate == "from species"),
              " species of our data set, and in addition allowed us to get a proxy of Ne for ",
              sum(data1$Ne_estimate == "from genus")," species",sep=""))
  
  
  print(paste("Among the ",nrow(data1)," species analyzed",sep=""))
  
  
  print(paste("C. elegans (S^{hx} =",
              round(data1["Caenorhabditis_elegans",]$S_POCs,1),") and C. nigoni (S^{hx} =",
              round(data1["Caenorhabditis_nigoni",]$S_POCs,1),")",sep=""))
  
  
  print(paste("Dipters also show relatively strong values of $S^{hx}$ (mean=",
              round(mean(data1[data1$clade_group == "Diptera",]$S_POCs),2),
              "±",round(sd(data1[data1$clade_group == "Diptera",]$S_POCs),2),
              " sd)",sep=""))
  
  
  print(paste("lepidopters (mean S^{hx}=",
              round(mean(data1[data1$clade_group == "Lepidoptera",]$S_POCs),2),
              "±",round(sd(data1[data1$clade_group == "Lepidoptera",]$S_POCs),2),
              " sd)",sep=""))
  
  print(paste("In vertebrates, signals of translational selection are weak (mean=",
              round(mean(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$S_POCs),2),
              "±",round(sd(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$S_POCs),2),
              " sd)",sep=""))
  
  ttest = t.test(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$S_POCs)
  
  print(paste("(Student's t-Test, p-value<",format(ttest$p.value,scientific = T),sep=""))
  
  
  ##  In species subject to translational selection, the tRNA pool evolves in response to changes in neutral substitution patterns
  print("In species subject to translational selection, the tRNA pool evolves in response to changes in neutral substitution patterns")
  
  
  print(paste("range of variation in GC-content (GCi ranging from ",round(min(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$gci),2),
              " to ",round(max(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$gci),2),",",sep=""))
  
  
  print(paste("that strongly correlates with their average GC3 (from ",round(min(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$gc3),2),
              " to ",round(max(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$gc3),2),",",sep=""))
  
  
  print(paste("we focused our analyses on the ",
              nrow(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]),
              " Diptera and Lepidoptera",sep=""))
  
  data16 = read.delim("data/data16_supp.tab",comment.char = "#")
  data16 = data16[data16$species %in% data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$species,]
  nna_nng_code = data16[substr(data16$codon,3,3) %in% c("A","G") & !data16$amino_acid %in% c('Ter','Trp','Met','Ile'),]
  
  print(paste("we observed that the decoding of ",sum(tapply(nna_nng_code$nb_tRNA_copies == 0,nna_nng_code$aa_name_scu,sum) ==0) , " NNA/NNG synonymous codon pairs",sep=""))
  
  
  ##  Variation in optimality between wobble and Watson-Crick pairing
  print("Variation in optimality between wobble and Watson-Crick pairing")
  
  
  duc_nna_nng = nna_nng_code[order(nna_nng_code$rank),]
  duc_nna_nng = duc_nna_nng[paste(duc_nna_nng$species,duc_nna_nng$aa_name_scu) %in% 
                              paste(duc_nna_nng[duc_nna_nng$nb_tRNA_copies == 0,]$species,duc_nna_nng[duc_nna_nng$nb_tRNA_copies == 0,]$aa_name_scu) ,]
  duc_nna_nng = duc_nna_nng[!duplicated(paste(duc_nna_nng$species,duc_nna_nng$aa_name_scu)),]
  
  nnt_nnc_code = code[substr(code$codon,3,3) %in% c("T","C"),]
  print(paste("NNT/NNC synonymous codon pairs (N=" , nrow(nnt_nnc_code)/2 , "),",sep=""))
  
  data_nnc_nnu = data16[substr(data16$codon,3,3) %in% c("C","T"),]
  
  print(paste("for ",
              round(nrow(data_nnc_nnu[data_nnc_nnu$nb_tRNA_copies == 0,])/(nrow(data_nnc_nnu)/2),2)*100,"% of the ",nrow(data_nnc_nnu)/2,
              " NNT/NNC synonymous codons pairs analyzed (" , nrow(nnt_nnc_code)/2 , " pairs $\times$ ",
              length(unique(data_nnc_nnu$species))," species) ",sep=""))
  
  
  print(paste("the ",length(c("Phe", "Tyr", "Cys", "His", "Asn", "Asp","Ser_2")),
              " pairs corresponding to the amino acids with duet codons (Phe, Tyr, Cys, His, Asn, Asp and the AGT/AGC 'duet' of Ser), and the ",
              length(c("Ile","Ser", "Leu", "Pro", "Arg", "Thr", "Val", "Ala", "Gly"))," pairs from amino acids with triplet (Ile) or quartet codons (Ser, Leu, Pro, Arg, Thr, Val, Ala, Gly).",sep=""))
  
  
  duet_nnc_nnu = data_nnc_nnu[data_nnc_nnu$aa_name_scu %in% c("Phe", "Tyr", "Cys", "His", "Asn", "Asp","Ser_2"),]
  duc_nnc_nnu = duet_nnc_nnu[order(duet_nnc_nnu$rank),]
  duc_nnc_nnu = duc_nnc_nnu[paste(duc_nnc_nnu$species,duc_nnc_nnu$aa_name_scu) %in% 
                              paste(duc_nnc_nnu[duc_nnc_nnu$nb_tRNA_copies == 0,]$species,duc_nnc_nnu[duc_nnc_nnu$nb_tRNA_copies == 0,]$aa_name_scu) ,]
  duc_nnc_nnu = duc_nnc_nnu[!duplicated(paste(duc_nnc_nnu$species,duc_nnc_nnu$aa_name_scu)),]
  
  print(paste("For NNC/NNT duets, when a single tRNA is present (",
              round(sum(duet_nnc_nnu$nb_tRNA_copies == 0) / 
                      (nrow(duet_nnc_nnu)/2),2)*100,
              "% of cases), it is always the GNN-tRNA, and in ",
              round(table(substr(duc_nnc_nnu$codon,3,3))['C']/(nrow(duc_nnc_nnu)),2)*100,
              "% of cases it is the NNC codon, decoded through Watson-Crick pairing that is preferred.",sep=""))
  
  table(substr(duet_nnc_nnu[duet_nnc_nnu$nb_tRNA_copies == 0,]$codon,3,3)) # 'it is always the GNN-tRNA'
  
  other_pairs = data_nnc_nnu[!data_nnc_nnu$aa_name_scu %in% c("Phe", "Tyr", "Cys", "His", "Asn", "Asp","Ser_2"),]
  duc_other = other_pairs[order(other_pairs$rank),]
  duc_other = duc_other[paste(duc_other$species,duc_other$aa_name_scu) %in% 
                          paste(duc_other[duc_other$nb_tRNA_copies == 0,]$species,duc_other[duc_other$nb_tRNA_copies == 0,]$aa_name_scu) ,]
  duc_other = duc_other[!duplicated(paste(duc_other$species,duc_other$aa_name_scu)),]
  
  print(paste("For 8 of the 9 other pairs, when a single tRNA is present (",
              round(sum(other_pairs$nb_tRNA_copies == 0) / 
                      (nrow(other_pairs)/2),2)*100,"% of cases), it is always the ANN-tRNA, the only exception being Gly (GNN-tRNA).",sep=""))
  
  table(substr(other_pairs[other_pairs$nb_tRNA_copies == 0,]$codon,1,3)) # 'it is always the ANN-tRNA'
  
  print(paste("For Gly, the GGT codon, decoded via wobble pairing, is preferred to the GGC codon in ",
              round(table(duc_other[duc_other$aa_name_scu == "Gly",]$codon)["GGT"] / nrow(duc_other[duc_other$aa_name_scu == "Gly",]),2) * 100,"% of species.",sep=""))
  
  
  print(paste("For the other pairs (decoded by ANN-tRNA) there is more variability: the NNC codon (wobble pairing) is preferred in ",
              round(table(substr(duc_other[duc_other$aa_name_scu != "Gly",]$codon,3,3))["C"] / nrow(duc_other[duc_other$aa_name_scu != "Gly",]),2)*100,
              "% of species, whereas the NNT codon (watson-crick pairing) is preferred in the others.",sep=""))
  
  
  
  ##  Discussion
  print("Discussion")
  ##  Predicting translationally optimal codons
  print("Predicting translationally optimal codons")

  
  print(paste("higher than ours (respectively ",round(data1[data1$species == "Homo_sapiens",]$S_POCs,2)," and ",
              round(data1[data1$species == "Mus_musculus",]$S_POCs,2),")",sep=""))
  
  
  ##  Variation in the intensity of selection in favor of translationally optimal codons across metazoans
  print("Variation in the intensity of selection in favor of translationally optimal codons across metazoans")

  
  print(paste("Across the ",nrow(data1)," species, the highest values of $S$ are observed in \textit{Caenorhabditis} nematodes ($S=",round(data1[data1$species == "Caenorhabditis_elegans",]$S_POCs,2),"$ ",sep=""))
  
  
  print(paste(" and $S=",round(data1[data1$species == "Caenorhabditis_nigoni",]$S_POCs,2),"$ in \textit{C. nigoni}",sep=""))
  
  
  print(paste("We also found a clear signal of translational selection in diptera (mean $S=",
              round(mean(data1[data1$clade_group == "Diptera",]$S_POCs),2),"$, N=",nrow(data1[data1$clade_group == "Diptera",])," species)",sep=""))
  
  
  print(paste("lesser extent in lepidoptera (mean $S=",
              round(mean(data1[data1$clade_group == "Lepidoptera",]$S_POCs),2),"$, N=",nrow(data1[data1$clade_group == "Lepidoptera",])," species)",sep=""))
  
  
  print(paste("The weakness of translational selection in vertebrates (mean $S=",
              round(mean(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$S_POCs),2),
              "$, N=",nrow(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),])," species)",sep=""))
  
  
  
  print(paste("our dataset included ",nrow(data1[!data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei","Diptera") &
                                                   !grepl("Caenorhabditis",data1$species),])," invertebrate species",sep=""))
  
  
  
  ##  Materials & Methods
  print("Materials & Methods")
  ##  tRNA gene annotation
  print("tRNA gene annotation")
  
  data1 = read.delim("data/data1_supp.tab",comment.char = "#")
  print(paste("annotation file (N=",sum(!data1$tRNA_GFF)," species). For species in which tRNA annotations were not available (N=",sum(data1$tRNA_GFF)," species)",sep=""))
  
  
  
  # Supplementary Text
  print("Supplementary Text")
  
  ##  Highly constrained amino acids are enriched in optimals codons
  print("Highly constrained amino acids are enriched in optimals codons")
  
  data6 = read.delim("data/data6_supp.tab",comment.char = "#")
  dt_analysis = data6[data6$species == "Caenorhabditis_elegans",]
  
  print(paste("most constrained sites within proteins (from ",
              round(mean(dt_analysis[dt_analysis$set == "POCs" & dt_analysis$categorie == "Unconstrained",]$freq,na.rm=T)*100,1),"% on average to average ",
              round(mean(dt_analysis[dt_analysis$set == "POCs" & dt_analysis$categorie == "Highly constrained",]$freq,na.rm=T)*100,1),"%),",sep=""))
  
  
  print(paste("Overall, Delta POCcons is positive in ",round(sum(data1$constraint_overused_POCs>0)/nrow(data1)*100,),
              "% of species ",sep=""))
  
  
}
