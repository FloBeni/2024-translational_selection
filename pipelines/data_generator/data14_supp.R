# Generate Data 14
options(stringsAsFactors = F, scipen = 999)
library(stringr)
library(stringi)
library(tidyr)
library(dplyr)

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

subst_or_snp = "substitutions"
count_file = "subst"
ylabel = "Substitution rate"
proportion = 0.02

path_data = "/home/fbenitiere/data/"

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])

## Gene to study
GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species
species = "Drosophila_melanogaster"
print(species)
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid
path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")

codon_usage = read.delim( paste(path,"/codon_usage_gene_fpkm.txt.gz",sep="") )

codon_usage$length = rowSums(codon_usage[ , 3:66]) * 3
codon_usage$intern_stop_codon = rowSums(codon_usage[,stop_codon]) - grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon)

codon_usage = codon_usage[codon_usage$intern_stop_codon == 0 & codon_usage$start_codon == "ATG" & codon_usage$length_cds %% 3 == 0,]

if (quantile(grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),0.75) != 0){  # if annotated seq have a stop codon for the majority then remove those that dont
  codon_usage = codon_usage[grepl(paste(stop_codon,collapse = "|"),codon_usage$end_codon),] } else { print(species)}

codon_usage = codon_usage[order(codon_usage$length_cds,decreasing = T),]
codon_usage = codon_usage[!duplicated(codon_usage$gene_id),]
codon_usage = codon_usage[!is.na(codon_usage$median_fpkm) ,]
codon_usage = codon_usage[codon_usage$median_fpkm != 0 ,]
####

#### Gene to filter due to alignment problem
chrList = c(chrX="NC_004354.4",
            chr2L="NT_033779.5",chr2R="NT_033778.4",
            chr3L="NT_037436.4",chr3R="NT_033777.3",
            chr4="NC_004353.4", chrY="NC_024512.1")

treshold_list =  c("chr2L"="down_20000000", "chr2R"="up_7000000","chr4"="up_0",
                   "chr3L"="down_22000000","chr3R"="up_5000000","chrX"="down_20000000")
treshold_data = data.frame(treshold = sapply(treshold_list,function(x) as.numeric(str_split(x,'_')[[1]][2])))
treshold_data$chromosome = rownames(treshold_data)
treshold_data$orientation =  sapply(treshold_list,function(x) str_split(x,'_')[[1]][1])
rownames(treshold_data) = chrList[treshold_data$chromosome ] 


gene_data = read.delim("/home/fbenitiere/data/Projet-SplicedVariants/Annotations/Drosophila_melanogaster/formatted_data/gene.tab")
chrList = c("NC_004354.4"="chrX",
            "NT_033779.5" = "chr2L","NT_033778.4" = "chr2R",
            "NT_037436.4" = "chr3L","NT_033777.3" = "chr3R",
            "NC_004353.4" = "chr4", "NC_024512.1" = "chrY")
gene_data$chromosome = chrList[ gene_data$seq_id ]
gene_data$gene_id = sapply(gene_data$attributes ,function(x) str_replace(str_split(x,";")[[1]][1],'ID=',''))
gene_data = gene_data[ gene_data$chromosome %in% treshold_data$chromosome ,]

gene_data$filter_treshold = apply(gene_data,1,function(x){
  treshold = treshold_data[x["seq_id"],] 
  if (treshold$orientation == "down"){
    return(as.numeric(x["start"]) <= treshold$treshold & as.numeric(x["end"]) <= treshold$treshold   )
  }else if (treshold$orientation == "up"){
    return(as.numeric(x["start"]) >= treshold$treshold &as.numeric(x["end"]) >= treshold$treshold    )}
})

gene_data = gene_data[gene_data$filter_treshold,]

codon_usage = codon_usage[ codon_usage$gene_id %in% gene_data$gene_id ,]
#####

collect_subst_rate <- function(xaxis,intervalle_list,list_aa,list_codon){
  dt = data.frame(fpkm = tapply(xaxis, intervalle_list, median))
  for (region in c("codon","intron")){
    if(region == "codon"){
      df_synonymous = substitution_codon[substitution_codon$protein_id %in% count_codon_trinucl$protein_id,]
    } else if (region == "intron"){
      df_synonymous = substitution_trinucl[substitution_trinucl$protein_id %in% count_codon_trinucl$protein_id,]
    }
    df_synonymous = df_synonymous[ df_synonymous$aa_ancestral == df_synonymous$aa_subst &
                                     df_synonymous$aa_ancestral %in% list_aa &
                                     df_synonymous$ancestral %in% unlist(list_codon) & 
                                     df_synonymous$substituted %in% unlist(list_codon),]
    
    df_synonymous$codon = paste(
      names(list_codon)[sapply(df_synonymous$ancestral,function(x) grep(x,list_codon))],'->',
      names(list_codon)[sapply(df_synonymous$substituted,function(x) grep(x,list_codon))]
      ,sep="")
    
    for ( categorie in unique( df_synonymous$codon )){
      count_codon_trinucl[,paste(categorie,region,sep="_")] = table( df_synonymous[ df_synonymous$codon == categorie ,]$protein_id )[count_codon_trinucl$protein_id]
      if(region == "codon"){
        if (length(list_codon[str_split(categorie,"->")[[1]][1]][[1]]) != 1){
          count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] = rowSums( count_codon_trinucl[,list_codon[str_split(categorie,"->")[[1]][1]][[1]]] )
        } else {
          count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] = count_codon_trinucl[,list_codon[str_split(categorie,"->")[[1]][1]][[1]]]
        }
      } else if (region == "intron"){
        if (length(list_codon[str_split(categorie,"->")[[1]][1]][[1]]) != 1){
          count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] =  rowSums( count_codon_trinucl[,paste(list_codon[str_split(categorie,"->")[[1]][1]][[1]],"_intronic",sep="")] )
        } else {
          count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] =  ( count_codon_trinucl[,paste(list_codon[str_split(categorie,"->")[[1]][1]][[1]],"_intronic",sep="")] )
        }
      }
    }
    for ( categorie in unique( df_synonymous$codon )){
      name_col = str_replace(paste("density",categorie,region,sep="_"),"->","_to_")
      
      dt[,paste("sum_cible",name_col,sep="_")] = tapply( count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] , intervalle_list,function(x)  sum(x,na.rm=T))
      dt[,paste("sum_subst",name_col,sep="_")] = tapply( count_codon_trinucl[,paste(categorie,region,sep="_")] , intervalle_list,function(x)  sum(x,na.rm=T))
      
      dt[,name_col] = tapply( count_codon_trinucl[,paste(categorie,region,sep="_")] / count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] , intervalle_list,function(x)  mean(x,na.rm=T))
      
      
      dt_inter_err = data.frame()
      for (i in 1:100){ # bootstrapping
        avg_list = tapply( count_codon_trinucl[,paste(categorie,region,sep="_")] / count_codon_trinucl[,paste(str_split(categorie,"->")[[1]][1],region,sep="_")] , intervalle_list,function(x)  mean(sample(x,replace=T),na.rm=T))
        dt_inter_err = rbind(dt_inter_err,data.frame(
          average = avg_list,
          intervalle = names(avg_list)
        ))
      }
      
      
      dt[,paste("confint_low_",name_col,sep="")] = tapply(dt_inter_err$average,dt_inter_err$intervalle,function(x) quantile(x,c(0.025)))[rownames(dt)]
      dt[,paste("confint_high_",name_col,sep="")] = tapply(dt_inter_err$average,dt_inter_err$intervalle,function(x) quantile(x,c(0.975))) [rownames(dt)]
    }
  }
  return(dt)
}





list_file = list.files(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",sep=""),
                       pattern=paste("count",count_file,"_codon",sep=""))
analyzable_codon_count = data.frame()
for (file in list_file){
  data_per_chr = read.delim(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/" , file , sep=""))
  analyzable_codon_count = rbind(analyzable_codon_count,data_per_chr)
}

list_file = list.files(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",sep=""),
                       pattern=paste("count",count_file,"_trinucl",sep=""))
analyzable_trinucl_count = data.frame()
for (file in list_file){
  data_per_chr = read.delim(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",file,sep=""))
  analyzable_trinucl_count = rbind(analyzable_trinucl_count,data_per_chr)
}

analyzable_trinucl_count = analyzable_trinucl_count %>%
  mutate(protein_id = strsplit(as.character(protein_id), ";")) %>%
  unnest(protein_id) 


analyzable_trinucl_count = aggregate(analyzable_trinucl_count[,8:71], by=list(protein_id=analyzable_trinucl_count$protein_id), FUN=sum)
count_codon_trinucl = merge(x=analyzable_codon_count,y=analyzable_trinucl_count,by.x="protein_id",by.y = "protein_id",all.x=T,all.y=T,suffixes=c("","_intronic"))

count_codon_trinucl = count_codon_trinucl[count_codon_trinucl$protein_id %in% codon_usage$protein_id,]


list_file = list.files(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",sep="") , 
                       pattern=paste(count_file,"_cds",sep=""))

substitution_codon = data.frame()
for (file in list_file){
  data_per_chr = read.delim(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",file,sep=""))
  substitution_codon = rbind(substitution_codon,data_per_chr)
}
substitution_codon = substitution_codon[substitution_codon$protein_id %in% codon_usage$protein_id,]

substitution_codon$ancestral = substitution_codon$ancestral_codon
if (subst_or_snp == "snps" ){
  substitution_codon$substituted = substitution_codon$derived_codon
} else {
  substitution_codon$substituted = substitution_codon$substituted_codon
}

substitution_codon$aa_ancestral = code[ substitution_codon$ancestral,]$aa_name
substitution_codon$aa_subst = code[ substitution_codon$substituted,]$aa_name



list_file = list.files(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",sep=""),
                       pattern=paste(count_file,"_intron",sep=""))
substitution_trinucl = data.frame()
for (file in list_file){
  data_per_chr = read.delim(paste(path_data,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/",subst_or_snp,"/",file,sep=""))
  substitution_trinucl = rbind(substitution_trinucl,data_per_chr)
}
substitution_trinucl = substitution_trinucl %>%
  mutate(protein_id = strsplit(as.character(protein_id), ";")) %>%
  unnest(protein_id) 

substitution_trinucl = substitution_trinucl[substitution_trinucl$protein_id %in% codon_usage$protein_id,]

substitution_trinucl$ancestral = substitution_trinucl$ancestral_triplet
if (subst_or_snp == "snps" ){
  substitution_trinucl$substituted = substitution_trinucl$derived_triplet
} else {
  substitution_trinucl$substituted = substitution_trinucl$substituted_triplet
}

substitution_trinucl$site_polymorph = apply(substitution_trinucl,1,function(x) {
  sum(str_split(x["ancestral"],"")[[1]] != str_split(x["substituted"],"")[[1]])
})
substitution_trinucl$aa_ancestral = code[ substitution_trinucl$ancestral,]$aa_name
substitution_trinucl$aa_subst = code[ substitution_trinucl$substituted,]$aa_name



print(paste("Number of analyzable codons = ",sum(count_codon_trinucl[,c(4:67)])))
print(paste("Number of codons containing substitutions = ",nrow(substitution_codon)))
print(paste("Number of analyzable triplets = ",sum(count_codon_trinucl[,c(68:131)],na.rm = T)))
print(paste("Number of triplets containing substitutions = ",nrow(substitution_trinucl)))


rownames(codon_usage) = codon_usage$gene_id
count_codon_trinucl$median_fpkm =  codon_usage[count_codon_trinucl$gene_id,]$median_fpkm


species = "Drosophila_melanogaster"
print(species)
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid

path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")

tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
rownames(tRNA_optimal) = tRNA_optimal$codon
tRNA_optimal = tRNA_optimal[tRNA_optimal$aa_name != "Ter",]

subset_selected = tRNA_optimal[tRNA_optimal$POC2 | tRNA_optimal$POC1,]
list_aa = subset_selected$aa_name
list_POC = subset_selected$codon
list_nonPOC = tRNA_optimal[tRNA_optimal$aa_name %in% list_aa & !tRNA_optimal$codon %in% list_POC,]$codon

print(list_POC)
print(list_nonPOC)

list_codon = list(
  nonoptimal=list_nonPOC,
  optimal=list_POC
)

xaxis = count_codon_trinucl$median_fpkm
quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
intervalle_list = cut(xaxis, quantile,include.lowest = T,include.higher=T)

data14 = collect_subst_rate(xaxis,intervalle_list,list_aa,list_codon)

write.table(data14,"data/data14_supp.tab",quote=F,row.names = F,sep="\t")

