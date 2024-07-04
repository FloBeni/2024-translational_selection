# Generate Data 12

code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab",comment.char = "#")
rownames(GTDrift_list_species) = GTDrift_list_species$species


species_list = c( "Bactrocera_oleae","Ceratitis_capitata" , "Hermetia_illucens","Aedes_aegypti"  )

data12 = data.frame()
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
  
  
  xaxis = codon_usage$median_fpkm 
  proportion = 2/100
  quantile = unique( quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T ))
  intervalle_FPKM = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  FPKM_bins = tapply(xaxis, intervalle_FPKM, median)
  
  for (aa_name in unique(code$aa_name)){
    POC_codon = rowSums(codon_usage[ code[code$aa_name == aa_name,]$codon ],na.rm = T)
    for (codon in code[code$aa_name == aa_name,]$codon ){
      POC_obs = rowSums(codon_usage[ codon ],na.rm = T)
      
      data12 = rbind(data12,data.frame(
        species,
        aa_name,
        trinucl = substr(codon,3,3),
        trinucl2 = substr(codon,1,2),
        codon = codon,
        nb_genes = as.numeric(table(intervalle_FPKM)),
        rscu =  tapply( POC_obs / POC_codon*length(code[code$aa_name=="Val",]$codon)  , intervalle_FPKM , function(x) mean(x,na.rm=T)),
        fpkm = FPKM_bins
      ))
    }
  }
}

write.table(data12 , "data/data12_supp.tab",quote=F,row.names = F,sep="\t")
