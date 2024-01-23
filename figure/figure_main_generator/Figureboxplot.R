
tRNA_abundance = read.delim("/home/fbenitiere/data/papers/2024-translational_selection/data/tRNA_abundance.tab")

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon
code$nb_syn = table(code$aa_name)[code$aa_name]
code$nb_syn_scu = table(paste(code$aa_name, substr(code$codon,1,2),sep="_"))[paste(code$aa_name, substr(code$codon,1,2),sep="_")]
code$anticodon = sapply(code$codon,function(x) chartr("TUACG","AATGC",stri_reverse(x))  )
code$aa_name_scu = code$aa_name
code[code$nb_syn == 6,]$aa_name_scu =  paste(code[code$nb_syn == 6,]$aa_name ,code[code$nb_syn == 6,]$nb_syn_scu  ,sep="_")

code = code[!code$aa_name %in% c("Ter","Met","Trp"),]

rownames(code) = code$anticodon

tRNA_abundance_data = data.frame()
for (anticodon in code$anticodon){
  dt = data.frame(abundance = tRNA_abundance[,anticodon])
  dt$species = rownames(tRNA_abundance)
  dt$amino_acid = code[anticodon,]$aa_name
  dt$anticodon = anticodon
  dt$codon = code[anticodon,]$codon
  tRNA_abundance_data = rbind(tRNA_abundance_data,dt)
}

tRNA_abundance_data$color = sapply(tRNA_abundance_data$codon,function(x)substr(x,3,3))

tRNA_abundance_data$amino_acid = factor(tRNA_abundance_data$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa_name))
tRNA_abundance_data$codon = str_replace_all(tRNA_abundance_data$codon,'T','U')

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
tRNA_abundance_data$codon = factor(tRNA_abundance_data$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
library(ggplot2)

set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))

set_color = c(A="#A6CEE3",T="#1F78B4",C="#B2DF8A",G="#33A02C")

pbxplot = ggplot(tRNA_abundance_data,aes(x=codon,y=abundance)) + geom_boxplot(aes(fill=color)) +  coord_cartesian(ylim=c(0,50)) +
  scale_fill_manual("",values = set_color) + facet_wrap(~amino_acid,scales = "free")+
  theme_bw() + theme(
    title =  element_text(size=36, family="economica"),
    
    legend.text =  element_text(size=40, family="economica"),
    strip.text = element_text(size=30, family="economica",face="bold"),
    legend.spacing.x = unit(1, 'cm'),
    legend.position="top",
    legend.box.spacing = unit(2, "cm"),
    axis.text.x =  element_text( size=15),
    axis.text.y =  element_text( size=15),
    plot.title = element_text(hjust = 0.5,margin = margin(0,0,20,0))
  ) + theme(legend.position='none')
pbxplot

jpeg(paste(path_pannel,"pbxplotzoom.jpg",sep=""), width =2000/1,  1200/1,res=100/1)
print(pbxplot)
dev.off()
