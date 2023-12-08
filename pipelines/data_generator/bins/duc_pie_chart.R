# Generate Data 1
library(stringi)

path_data = "/home/fbenitiere/data/"

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])
code$nb_syn = table(code$aa_name)[code$aa_name]


GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

data1 = read.delim("data/data1.tab")
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1 = data1[data1$clade_group == "Mecopterida",]

data = data.frame()
for (species in data1$species){
  print(species)
  genome_assembly = GTDrift_list_species[species,]$assembly_accession
  taxID = GTDrift_list_species[species,]$NCBI.taxid
  
  path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")
  
  codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.tab.gz",sep="") )
  
  codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
  codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)
  
  codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]
  
  if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
    codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}
  
  codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
  codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
  codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
  codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
  
  observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )
  
  
  
  
  tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
  rownames(tRNA_optimal) = tRNA_optimal$codon
  
  xaxis = codon_usage$median_fpkm 
  proportion = 5/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  data_optiplus = data.frame()
  for ( amino_acid in unique(tRNA_optimal$aa_name)){
    codon_used = rownames(tRNA_optimal[tRNA_optimal$aa_name == amino_acid,])
    amino_acid_count = rowSums(codon_usage[codon_used],na.rm = T)
    triplet_intronic = rowSums(codon_usage[paste(codon_used,'_intronic',sep = "")],na.rm = T)
    
    for ( codon in codon_used ){
      COA_obs =  unlist(codon_usage[codon])
      COA_neg_obs = amino_acid_count
      
      COA_obs_intronic =  unlist(codon_usage[paste(codon,'_intronic',sep = "")])
      COA_neg_obs_intronic = triplet_intronic
      
      data_optiplus = rbind(data_optiplus,data.frame(
        species,
        amino_acid,
        codon,
        expressed_overused_background = (round(tapply( COA_obs / COA_neg_obs , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) - 
                                           round( mean( (COA_obs / COA_neg_obs)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T) , 5)) - (
                                             round(tapply( COA_obs_intronic / COA_neg_obs_intronic   , intervalle_FPKM , function(x) mean(x,na.rm=T))[length(FPKM_bins)],5) -
                                               round( mean( (COA_obs_intronic / COA_neg_obs_intronic)[codon_usage$median_fpkm <= median(codon_usage$median_fpkm )] , na.rm=T),5)
                                           )
      ))
    }
  }
  data_optiplus = data_optiplus[order(data_optiplus$expressed_overused_background,decreasing = T),]
  data_optiplus = data_optiplus[order(data_optiplus$amino_acid,decreasing = F),]
  data_optiplus$rank = unlist(tapply(data_optiplus$expressed_overused_background,data_optiplus$amino_acid,function(x)rev(rank(x))))
  data = rbind(data,data_optiplus)
}





data4 = data.frame()
for ( codon in unique(code$codon)){
  total = nrow(data[data$codon == codon,])
  # dc = data.frame(table(data[data$codon == codon,"rank"]))
  dc = data.frame(table(data[data$codon == codon,"expressed_overused_background"] > 0))
  dc$Prop = dc$Freq / total
  dc$total = total
  dc$amino_acid = paste(code[codon,]$aa_name , " (",code[codon,]$nb_syn,")",sep="")
  dc$codon = codon
  data4 = rbind(data4,dc)
}


data4$rank = factor(data4$Var1,levels =  c(1:6))
data4$Var1 = as.character(data4$Var1)
data4[data4$rank == 1,"categorie"] = "Best DUC"
data4[data4$Var1 == "TRUE","categorie"] = "DUC"
data4[data4$Var1 == "FALSE","categorie"] = "not DUC"
data4$categorie = factor(data4$categorie,levels =  c("not DUC","DUC"))
dt_graph = data4



vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
dt_graph$codon = factor(dt_graph$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","T","A","G"),sep=""))) )

dt_graph$title = factor(paste(dt_graph$codon,sep=""),  
                        sapply(levels(dt_graph$codon),function(x) paste(x,sep="")) )

dt_graph[dt_graph$amino_acid == "Ter (3)",]$Prop = NA
dt_graph[dt_graph$amino_acid == "Met (1)",]$Prop = NA
dt_graph[dt_graph$amino_acid == "Trp (1)",]$Prop = NA

{
p3A = ggplot(dt_graph, aes(x = "" , y = Prop, fill = fct_inorder(categorie))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") + facet_wrap(~title+paste(amino_acid),nrow=4,dir="v")+
  theme_void() + theme(
    title =  element_text(size=36, family="economica"),
    legend.text =  element_text(color="black", size=40, family="economica",vjust = 1.5,margin = margin(t = 10)),
    strip.text = element_text(size=30, family="economica",face="bold"),
    legend.spacing.x = unit(1, 'cm'),
    legend.position="right",
    legend.key.size= unit(1, "cm"),
    legend.box.margin=margin(l=50),
    plot.title = element_text(hjust = 0.5,margin = margin(0,0,20,0))
  )  +
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))+
  # scale_fill_manual("",values = c("#33A02C","#B2DF8A","#FF7F00","#FDBF6F","#FB9A99","#E31A1C"))  + 
  scale_fill_manual("",values = c("#FB9A99","#B2DF8A","#FF7F00","#FDBF6F","#FB9A99","#E31A1C"))
  # ggtitle(paste("Mecopterida"," (",dt_graph$total[1],")",sep=""))
p3A


jpeg(paste(path_pannel,"p3_DUC.jpg",sep=""), width = 10000/1,  3600/1,res=300/1)
print(p3A)
dev.off()
}
