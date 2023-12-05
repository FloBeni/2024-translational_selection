# Generate Data 7

options( stringsAsFactors = F, scipen = 999 )
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <-
    lapply(sheets, function(X)
      readxl::read_excel(filename, sheet = X))
  if (!tibble)
    x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
library(RColorBrewer)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))


########## 
mysheets <- read_excel_allsheets("/home/fbenitiere/data/Projet-NeGA/translational_selection/behrens_et_al/mmc2.xlsx")
names(mysheets)
data = mysheets[["Unique tRNA DE analysis"]]
colnames(data) = data[3,]
data = data[4:372,]
colnames(data)  = str_replace_all(colnames(data) ," ","_")
data = data[,!duplicated(colnames(data))]
data$iPSC_rep_2 = as.numeric(data$iPSC_rep_2)
data$iPSC_rep_1 = as.numeric(data$iPSC_rep_1)
data$Gene_copy_number = as.numeric(data$Gene_copy_number)

data$prop_iPSC_rep_2 = data$iPSC_rep_2 / sum(data$iPSC_rep_2) 
data$prop_iPSC_rep_1 = data$iPSC_rep_1 / sum(data$iPSC_rep_1) 
data$prop_Gene_copy_number = data$Gene_copy_number / sum(data$Gene_copy_number) 

p1=ggplot(data,aes(x=Gene_copy_number,y=prop_iPSC_rep_2)) + geom_smooth(method='lm', formula= y~x,col=set_color[1],size=2) +
  geom_point(pch=21,fill=set_color[2],size=3) + 
  ylab("Proportion uniquely aligned reads")+ xlab("Gene copy number") + 
  geom_abline() + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="serif"),
    axis.title.y = element_text(color="black", size=25, family="serif"),
    axis.text.y =  element_text(color="black", size=20, family="serif"),
    axis.text.x =  element_text(color="black", size=20,hjust=1, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.text =  element_text(color="black", size=18, family="serif")
  ) +   ggtitle(paste("H .sap, Rho Spearman =",round(cor.test( data$Gene_copy_number, data$prop_iPSC_rep_2,method="spearman",exact=F)$estimate,2)))
p1
cor.test( data$Gene_copy_number, data$prop_iPSC_rep_2,method="spearman",exact=F)$estimate
summary(lm(data$Gene_copy_number~ data$prop_iPSC_rep_2))

data$amino_acid = sapply(data$Unique_tRNA,function(x) str_split(x,'-')[[1]][1])

df = data.frame(
  species = "Homo_sapiens",
  amino_acid = names(tapply(as.numeric(data$Gene_copy_number),data$amino_acid,sum)),
  gene_copies = tapply(as.numeric(data$Gene_copy_number),data$amino_acid,sum),
  abundance_rep_1 = tapply(as.numeric(data$iPSC_rep_1),data$amino_acid,sum),
  abundance_rep_2 = tapply(as.numeric(data$iPSC_rep_2),data$amino_acid,sum)
)


GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species

#### Compare with codon usage
species = "Homo_sapiens"
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid
path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")


code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])
code$nb_syn = table(code$aa_name)[code$aa_name]

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

tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
rownames(tRNA_optimal) = tRNA_optimal$codon
tRNA_optimal = tRNA_optimal[tRNA_optimal$aa_name != "Ter",]

observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )

aa_data = data.frame()
for (amino_acid in unique(code$aa_name)){
  codon_used = rownames(code[code$aa_name == amino_acid,])
  aa_data = rbind(aa_data,
                  data.frame(
                    amino_acid ,
                    letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                    tRNASE_copies= sum(tRNA_optimal[codon_used,]$nb_tRNA_copies,na.rm = T),
                    obs_codon = sum(observation[codon_used])
                  ))
}
aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]

rownames(aa_data) = aa_data$amino_acid

df$amino_acid[df$amino_acid == "iMet"] = "Met"

df$transcriptome_count = aa_data[df$amino_acid,]$obs_codon
df = df[!is.na(df$transcriptome_count),]

df$prop_gene_copies = df$gene_copies / sum(df$gene_copies)
df$prop_abundance_average = (df$abundance_rep_1 / sum(df$abundance_rep_1) + df$abundance_rep_2 / sum(df$abundance_rep_2)) /2
df$prop_transcriptome_count = df$transcriptome_count / sum(df$transcriptome_count)


data7 = df


##########

dt_droso = read.table("/home/fbenitiere/data/Projet-NeGA/translational_selection/behrens_et_al/bg3-normCounts.csv",sep=",",header=T)
dt_droso$Anticodon = str_replace(dt_droso$Gene,"Drosophila_melanogaster_tRNA-","")
dt_droso$amino_acid = sapply(dt_droso$Anticodon,function(x) str_split(x,'-')[[1]][1])

dt_droso$prop_iPSC_rep_2 = dt_droso$lib250_2_dt_bg3_2_42_on.unpaired_uniq / sum(dt_droso$lib250_2_dt_bg3_2_42_on.unpaired_uniq)

p2=ggplot(dt_droso,aes(x=size,y=prop_iPSC_rep_2)) + geom_smooth(method='lm', formula= y~x,col=set_color[3],size=2) +
  geom_point(pch=21,fill=set_color[4],size=3) + 
  ylab("Proportion uniquely aligned reads")+ xlab("Gene copy number") + 
  geom_abline() + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="serif"),
    axis.title.y = element_text(color="black", size=25, family="serif"),
    axis.text.y =  element_text(color="black", size=20, family="serif"),
    axis.text.x =  element_text(color="black", size=20,hjust=1, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.text =  element_text(color="black", size=18, family="serif")
  ) +   ggtitle(paste("D .mel, Rho Spearman =",round(cor.test( dt_droso$size, dt_droso$prop_iPSC_rep_2,method="spearman",exact=F)$estimate,2)))
p2

summary(lm(dt_droso$lib250_2_dt_bg3_2_42_on.unpaired_uniq~dt_droso$size))
cor.test( dt_droso$lib250_2_dt_bg3_2_42_on.unpaired_uniq, dt_droso$size,method="spearman",exact=F)



df = data.frame(
  species = "Drosophila_melanogaster",
  amino_acid = names(tapply(dt_droso$lib250_1_dt_bg3_1_42_on.unpaired_uniq,dt_droso$amino_acid,sum)),
  gene_copies = tapply(dt_droso$size,dt_droso$amino_acid,sum),
  abundance_rep_1 = tapply(dt_droso$lib250_1_dt_bg3_1_42_on.unpaired_uniq,dt_droso$amino_acid,sum),
  abundance_rep_2 =tapply(dt_droso$lib250_2_dt_bg3_2_42_on.unpaired_uniq,dt_droso$amino_acid,sum)
)



#### Compare with codon usage D.melano
species = "Drosophila_melanogaster"
genome_assembly = GTDrift_list_species[species,]$assembly_accession
taxID = GTDrift_list_species[species,]$NCBI.taxid
path = paste("data/per_species/",species,"_NCBI.taxid",taxID,"/",genome_assembly,sep="")


code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
stop_codon = rownames(code[code$aa_name == "Ter",])
code$nb_syn = table(code$aa_name)[code$aa_name]

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

tRNA_optimal = read.delim(paste(path,"/decoding_table.tab.gz",sep=""))
rownames(tRNA_optimal) = tRNA_optimal$codon
tRNA_optimal = tRNA_optimal[tRNA_optimal$aa_name != "Ter",]

observation = colSums( codon_usage[3:70] * codon_usage$median_fpkm , na.rm = T )

aa_data = data.frame()
for (amino_acid in unique(code$aa_name)){
  codon_used = rownames(code[code$aa_name == amino_acid,])
  aa_data = rbind(aa_data,
                  data.frame(
                    amino_acid ,
                    letter_aa = unique(code[code$aa_name == amino_acid,]$aa) ,
                    tRNASE_copies= sum(tRNA_optimal[codon_used,]$nb_tRNA_copies,na.rm = T),
                    obs_codon = sum(observation[codon_used])
                  ))
}
aa_data = aa_data[!grepl("Ter",aa_data$amino_acid) ,]

rownames(aa_data) = aa_data$amino_acid

df$amino_acid[df$amino_acid == "iMet"] = "Met"

df$transcriptome_count = aa_data[df$amino_acid,]$obs_codon
df = df[!is.na(df$transcriptome_count),]

df$prop_gene_copies = df$gene_copies / sum(df$gene_copies)
df$prop_abundance_average = (df$abundance_rep_1 / sum(df$abundance_rep_1) + df$abundance_rep_2 / sum(df$abundance_rep_2)) /2
df$prop_transcriptome_count = df$transcriptome_count / sum(df$transcriptome_count)


data7 = rbind(data7,df)



write.table(data7,"data/data7_supp.tab",quote=F,row.names = F,sep="\t")
