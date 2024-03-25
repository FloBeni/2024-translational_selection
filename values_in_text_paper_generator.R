source("figure/figure_main_generator/library_path.R")


data1 = read.delim("data/data1_supp.tab")
rownames(data1) = data1$species
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group


code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code = code[code$aa_name != "Ter",]

##  Non-adaptive processes are the primary drivers of codon usage variations among metazoans
print("Non-adaptive processes are the primary drivers of codon usage variations among metazoans")

print(paste("position between ",nrow(data1)," metazoan species",sep=""))


print(paste("gene expression data (N=",sum(data1$nb_genes_filtered<5000)," species with less than 5,000 genes",sep=""))


print(paste("CU across ",sum(data1$nb_genes_filtered>=5000)," metazoan species",sep=""))

data1 = data1[data1$nb_genes_filtered>=5000,]

print(paste(sum(data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"))," vertebrates, ",
            sum(data1$clade_group %in% c("Diptera","Lepidoptera","Coleoptera","Hymenoptera","Other Insects")),
            " insects and ",sum(data1$clade_group %in%c("Other Metazoans","Nematoda"))," other metazoan clades",sep=""))

print(paste("constrained and spans from ",round(min(data1$gc3),2)," to ",round(max(data1$gc3),2),sep=""))

print(paste(" dipterans (N=",sum(data1$clade_group %in% c("Diptera")),") exhibit",sep=""))

print(paste("GC3 variations (from ",round(min(data1[data1$clade_group %in% c("Diptera"),]$gc3),2)," to ",round(max(data1[data1$clade_group %in% c("Diptera"),]$gc3),2),")",sep=""))


print(paste("GC3 and GCi are highly correlated (rho=",round(data1["Homo_sapiens",]$rho_gc3_gci,2),")",sep=""))

print(paste("pronounced in Caenorhabditis elegans (rho=",round(data1["Caenorhabditis_elegans",]$rho_gc3_gci,2),",",sep=""))

dt_graph = data1
ylabel = "var_gc3"
xlabel = "var_gci"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 
dt_graph[,c(ylabel,xlabel)] = sqrt(dt_graph[,c(ylabel,xlabel)])

model_to_use = fitted_model(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

print(paste("were higly correlated (R^2=",model_to_use$r2,")",sep=""))


print(paste(" gene GCi and GC3 (SD_{GCi} > ",
            round(min(sqrt(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Hymenoptera"),]$var_gci)),3),
            ", SD_{GC3} > ",
            round(min(sqrt(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Hymenoptera"),]$var_gc3)),3),sep=""))



##  tRNA abundance matches proteome requirements
print("tRNA abundance matches proteome requirements")


tRNA_abundance = read.delim("/home/fbenitiere/2024-translational_selection/data/tRNA_abundance.tab")
tRNA_abundance = tRNA_abundance[,code$anticodon]
tRNA_abundance$clade_group = GTDrift_list_species[rownames(tRNA_abundance),]$clade_group


print(paste("ranging from an average of ",
            round(mean(rowSums(tRNA_abundance[tRNA_abundance$clade_group == 'Hymenoptera',code$anticodon])),0),
            " tRNA gene copies per genome in hymenopterans to ",
            format(round(mean(rowSums(tRNA_abundance[tRNA_abundance$clade_group == 'Teleostei',code$anticodon])),0),big.mark=",",scientific=F)
            ," copies in teleost fish",sep=""))




data7 = read.delim("data/data7_supp.tab")
dt_graph = data7[data7$species == "Drosophila_melanogaster",]
spearman_method_aa = cor.test( dt_graph$prop_abundance_average, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
print(paste("direct measures of tRNA abundance (rho = ",
            round(spearman_method_aa$estimate,2),sep=""))
dt_graph = data7[data7$species == "Homo_sapiens",]
spearman_method_aa_human = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
dt_graph = data7[data7$species == "Drosophila_melanogaster",]
spearman_method_aa_droso = cor.test( dt_graph$gene_copies, dt_graph$prop_transcriptome_count,method="spearman",exact=F)
print(paste(" with their tRNA gene copy numbers (rho=",round(spearman_method_aa_droso$estimate,2)," and ",round(spearman_method_aa_human$estimate,2)," respectively",sep=""))


print(paste(" tRNA gene copy number also correlates with the amino acid demand in Caenorhabditis elegans (rho=",
            round(data1['Caenorhabditis_elegans',]$rho_aa_fpkm,2),sep=""))


print(paste(" same analysis conducted across ",nrow(data1)," animal species ",sep=""))

print(paste("in ",round(sum(data1$pval_aa_fpkm<=0.05)/nrow(data1),2)*100,"% of the species ",sep=""))

print(paste("correlates with amino acid usage (N=",sum(data1$pval_aa_fpkm<=0.05)," species)",sep=""))

data1 = data1[data1$pval_aa_fpkm<=0.05,]
print(paste("also excluded ",sum(data1$nb_codon_not_decoded!=0)," species ",sep=""))

data1 = data1[data1$nb_codon_not_decoded == 0,]

##  Definition of putative-optimal codons based on tRNA abundance and wobble-pairing rules
print("Definition of putative-optimal codons based on tRNA abundance and wobble-pairing rules")

print(paste(nrow(code)," codons to their cognate tRNA",sep=""))

tRNA_abundance = tRNA_abundance[rownames(tRNA_abundance)%in%data1$species,]
print(paste("ranges from ",min(rowSums(tRNA_abundance[,code$anticodon] != 0)),
            " to ",max(rowSums(tRNA_abundance[,code$anticodon] != 0)),
            " per species (average=",round(mean(rowSums(tRNA_abundance[,code$anticodon] != 0)),0),
            ").",sep=""))

print(paste("This implies that ",61-max(rowSums(tRNA_abundance[,code$anticodon] != 0)),
            " to ",61-min(rowSums(tRNA_abundance[,code$anticodon] != 0)),sep=""))


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
print(paste("AAT accounts for ",round(AAT/(AAC+AAT)*100),"% of asparagine codons, highlighting the significance of wobble pairing",sep=""))

print(paste("There are ",length(unique(code[code$nb_syn>1 & code$aa_name != "Ter",]$aa_name))," amino acids that",sep=""))



data12 = read.delim("data/data12_supp.tab")
data12 = data12[data12$species %in% data1$species,]

list_species = c()
list_exaequo = c()
for(species in unique(data12$species)){
  aa_present = table(data12[data12$species == species & data12$nb_syn != 1 & data12$abundance != 0 ,]$amino_acid)
  aa_present = names(aa_present[aa_present >= 2 ])
  if (any(!(aa_present %in% unique(data12[data12$species == species & data12$amino_acid %in% aa_present & data12$POC1, ]$amino_acid)))){
    list_species=append(list_species,species)
    list_exaequo=append(list_exaequo,aa_present[!(aa_present %in% unique(data12[data12$species == species & data12$amino_acid %in% aa_present & data12$POC1, ]$amino_acid))])
  }
}

print(paste(length(list_exaequo)," ex æquo cases in total, in ",length(list_species)," species",sep=""))

print(paste("In the case of human we are able to identify POC1 for ",length(unique(data12[data12$species == "Homo_sapiens" & data12$POC1,]$amino_acid)),
            " amino acid and POC2 for ",length(unique(data12[data12$species == "Homo_sapiens" & data12$POC2,]$amino_acid))," codons.",sep=""))



print(paste("For the human genome, POC1 have been defined for ",length(unique(data12[data12$species == "Homo_sapiens" & data12$POC1,]$amino_acid)),
            " amino acids and POC2 for ",length(unique(data12[data12$species == "Homo_sapiens" & data12$POC2,]$amino_acid)),
            " amino acids",sep=""))


print(paste("In contrast, for Caenorhabditis elegans, POC1 and POC2 are defined for ",length(unique(data12[data12$species == "Caenorhabditis_elegans" & data12$POC1,]$amino_acid)),
            " and ",length(unique(data12[data12$species == "Caenorhabditis_elegans" & data12$POC2,]$amino_acid)),
            " amino acids",sep=""))


print(paste("On average among the ",nrow(data1),
            " species, POC1 are defined for ", round(mean(data1$nb_aa_POC1),1),
            " amino acids per species (ranging from ",min(data1$nb_aa_POC1)," to ",max(data1$nb_aa_POC1),
            ") and POC2 for ", round(mean(data1$nb_aa_POC2),1),
            " amino acids (ranging from ",min(data1$nb_aa_POC2)," to ",max(data1[data1$species!="Tyto alba",]$nb_aa_POC2)-1,
            ", except Tyto alba with ",data1["Tyto_alba",]$nb_aa_POC2," including Ile)",sep=""))

##  Highly expressed genes are enriched in optimal codons
print("Highly expressed genes are enriched in optimal codons")

data5 = read.delim("data/data5_supp.tab")

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


print(paste("which represent ",round(sum(codon_usage$median_fpkm <= 50)/nrow(codon_usage),2)*100,"% of genes in \textit{C. elegans}",sep=""))


print(paste("the studied species (N=",nrow(data1)," species)",sep=""))

summ_lm = summary(lm(data1$expressed_overused_background_POC1~data1$expressed_overused_background_POC2))
print(paste("are strongly correlated (R2=",round(summ_lm$r.squared,2)*100,"% , p-value<",
            format(summ_lm$coefficients[2,4],scientific = T),sep=""))

print(paste(" For ",sum(data1$expressed_overused_background_POC1 >=0,na.rm = T )," species (",
            round(sum(data1$expressed_overused_background_POC1 >=0,na.rm = T )/nrow(data1),2)*100,"%),",sep=""))


print(paste("are around +",round(mean(data1[data1$clade_group == "Diptera",]$expressed_overused_background_POCs)),
            "% in Diptera compared to +",
            round(mean(data1[data1$clade_group %in% c("Mammalia","Aves","Other Tetrapods","Teleostei"),]$expressed_overused_background_POCs)),
            "% in vertebrates",sep=""))



##  Highly constrained amino acids are enriched in optimals codons
print("Highly constrained amino acids are enriched in optimals codons")

data6 = read.delim("data/data6_supp.tab")
dt_analysis = data6[data6$species == "Caenorhabditis_elegans",]

print(paste("most constrained sites within proteins (from ",
            round(mean(dt_analysis[dt_analysis$type_aa == "POCs" & dt_analysis$categorie == "Unconstrained",]$freq,na.rm=T)*100,1),"% on average to average ",
            round(mean(dt_analysis[dt_analysis$type_aa == "POCs" & dt_analysis$categorie == "Highly constrained",]$freq,na.rm=T)*100,1),"%),",sep=""))


print(paste("Overall, Delta POCcons is positive in ",round(sum(data1$constraint_overused_POCs>0)/nrow(data1)*100,),
            "% of species ",sep=""))


##  Selection favors optimal codons in highly expressed genes of Drosophila melanogaster
print("Selection favors optimal codons in highly expressed genes of Drosophila melanogaster")


data_11 = read.delim("data/data11_supp.tab")
print(paste(format(sum(data_11$sum_subst_density_nonoptimal_to_optimal_codon),big.mark=",",scientific=T),
            " nPO$>$PO synonymous SNPs and ",
            format(sum(data_11$sum_subst_density_optimal_to_nonoptimal_codon),big.mark=",",scientific=T)," PO$>$nPO synonymous SNPs",sep=""))


data_10 = read.delim("data/data10_supp.tab")
print(paste(format(sum(data_10$sum_subst_density_nonoptimal_to_optimal_codon),big.mark=",",scientific=T),
            " nPO$>$PO synonymous substitutions and ",
            format(sum(data_10$sum_subst_density_optimal_to_nonoptimal_codon),big.mark=",",scientific=T)," PO$>$nPO synonymous substitutions",sep=""))

print(paste("In introns, we observed ",
            format(sum(data_11$sum_subst_density_nonoptimal_to_optimal_intron),big.mark=",",scientific=T)," nPO$>$PO SNPs and ",
            format(sum(data_11$sum_subst_density_optimal_to_nonoptimal_intron),big.mark=",",scientific=T)," PO$>$nPO SNPs, ",
            format(sum(data_10$sum_subst_density_nonoptimal_to_optimal_intron),big.mark=",",scientific=T)," nPO$>$PO substitutions and ",
            format(sum(data_10$sum_subst_density_optimal_to_nonoptimal_intron),big.mark=",",scientific=T)," PO$>$nPO substitutions. ",sep=""))


##  Weak relationship between the strength of translational selection and the effective population size
print("Weak relationship between the strength of translational selection and the effective population size")

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


print(paste("This list included ",sum(lynch_dt$species %in% data1$species),
            " species of our data set, and in addition allowed us to get a proxy of Ne for ",
            sum(!is.na(data1$Ne))," species",sep=""))


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


##  The tRNA pool evolves in response to changes in genomic substitution patterns
print("The tRNA pool evolves in response to changes in genomic substitution patterns")

print(paste("Diptera and Lepidoptera, that show a strong translational selection compared to other metazoans (N=",
            nrow(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]),
            " species; we excluded",sep=""))

print(paste("range of variation in GC-content (GCi ranging from ",round(min(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$gci),2),
            " to ",round(max(data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$gci),2),",",sep=""))


data_codon = read.delim("data/data_codons.tab")
data_codon = data_codon[data_codon$species %in% data1[data1$clade_group %in% c("Diptera","Lepidoptera") & data1$species !="Eumeta_japonica",]$species,]
nna_nng_code = data_codon[substr(data_codon$codon,3,3) %in% c("A","G") & !data_codon$aa_name %in% c('Ter','Trp','Met','Ile'),]

tapply(nna_nng_code$nb_tRNA_copies == 0,nna_nng_code$codon,sum)

print(paste("NNA/NNG synonymous codon pairs (N=",length(table(substr(nna_nng_code$codon,1,2))) , " pairs)",sep=""))

duc_nna_nng = nna_nng_code[order(nna_nng_code$rank),]
duc_nna_nng = duc_nna_nng[paste(duc_nna_nng$species,duc_nna_nng$aa_name_scu) %in% 
                             paste(duc_nna_nng[duc_nna_nng$nb_tRNA_copies == 0,]$species,duc_nna_nng[duc_nna_nng$nb_tRNA_copies == 0,]$aa_name_scu) ,]
duc_nna_nng = duc_nna_nng[!duplicated(paste(duc_nna_nng$species,duc_nna_nng$aa_name_scu)),]
table(duc_nna_nng$codon)



nnt_nnc_code = code[substr(code$codon,3,3) %in% c("T","C"),]
print(paste("For NNT/NNC synonymous codon pairs (N=" , nrow(nnt_nnc_code)/2 , "),",sep=""))

data_nnc_nnu = data_codon[substr(data_codon$codon,3,3) %in% c("C","T"),]

print(paste("Indeed, among the ",nrow(data_nnc_nnu)/2," NNT/NNC synonymous codons pairs analyzed (" , nrow(nnt_nnc_code)/2 , " pairs $\times$ ",
            length(unique(data_nnc_nnu$species))," species) ",
            round(nrow(data_nnc_nnu[data_nnc_nnu$nb_tRNA_copies == 0,])/(nrow(data_nnc_nnu)/2),2)*100,"% are decoded by a ",sep=""))


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



print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
print(paste("",sep=""))
