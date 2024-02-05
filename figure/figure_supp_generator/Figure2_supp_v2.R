# Generate Supplementary Figure 2
source("figure/figure_supp_generator/library_path.R")

require(ade4)
require(adephylo)
require(ape)
require(ggtree)
require(treeio)
require(tidyr)
require(tibble)
require(ggplot2)
require(ggExtra)
require(aplot)
require(dplyr)
require(adegraphics)

tRNA <- read.table("data/tRNA_abundance.tab")
rownames(tRNA) = str_replace_all(rownames(tRNA) , "_"," ")
newick <- readLines('data/GTDrift_metazoa_phylogenetic_tree.nwk')

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$anticodon

# phylog <- newick2phylog(newick)
# table.phylog(tRNA, phylog, cnodes='', labels.row='', grid=F, cleaves=0, csize=1)
# tree <- read.tree(text = newick)
# p <- ggtree(tree)
arbrePhylotips = arbrePhylo
arbrePhylotips$tip.label <- str_replace_all(arbrePhylotips$tip.label,"_"," ")
edge_group <- str_replace_all(arbrePhylotips$tip.label,"_"," ")
edge_clade <- rep("branch",length(arbrePhylotips$edge[,2]))
for (group in unique(edge_group)){
  if (group %in% unlist(listNomSpecies)){
    edge_clade[arbrePhylotips$edge[,2] %in% grep(group,edge_group)] =
      names(listNomSpecies[unlist(lapply(listNomSpecies,function(x) group %in% x))])
  }
}

edge_clade_prev = edge_clade
list_inclusion =  list("Other Metazoans"=c("Diptera","Lepidoptera","Coleoptera","Other Insects","Other Tetrapods","Nematoda","Hymenoptera","Teleostei","Other Metazoans"),
                       "Other Tetrapods"=c("Aves","Mammalia","Other Tetrapods"),"Other Insects"=c("Diptera","Lepidoptera","Coleoptera","Other Insects","Hymenoptera"),
                       Nematoda="Nematoda",Teleostei="Teleostei",Hymenoptera="Hymenoptera",Aves="Aves",Mammalia="Mammalia","Diptera"="Diptera","Lepidoptera"="Lepidoptera","Coleoptera"="Coleoptera"
)



clade="Other Metazoans"

for (clade in names(list_inclusion)){
  edge_clade[ which.edge(arbrePhylotips,  arbrePhylotips$edge[,2][edge_clade_prev %in% unlist(list_inclusion[clade])] ) ] = clade
}
node_metadata = data.frame(node=arbrePhylotips$edge[,2],color=edge_clade)

node_metadata$color = factor(node_metadata$color, levels = names(Clade_color))

offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree") # to use flip deprecated

pA = ggtree(arbrePhylotips) %>% flip(264, 375)
pA <- pA %<+% node_metadata  + aes(color=color) + 
  scale_color_manual("Clade",values=Clade_color[unique(edge_clade)]) + theme(
    panel.background = element_rect(fill = "white", linetype = "dashed")
  )  + theme(legend.position = "none")

pA
# + geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,size=0)
# tidy the data
tRNA_long <- 
  tRNA %>% rownames_to_column(var='sp') %>%
  pivot_longer(cols=!sp,
               names_to = 'tRNA',
               values_to = "count") 

tRNA_long$clade_group = GTDrift_list_species[str_replace_all(tRNA_long$sp," ","_"),]$clade_group
tRNA_long$aa_name = code[tRNA_long$tRNA,]$aa_name
tRNA_long$amino_acid = code[tRNA_long$tRNA,]$aa_name
tRNA_long$anticodon = code[tRNA_long$tRNA,]$anticodon
tRNA_long$codon = code[tRNA_long$tRNA,]$codon
tRNA_long$amino_acid = factor(tRNA_long$aa_name,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa_name))
tRNA_long$anticodon = str_replace_all(tRNA_long$anticodon,'T','U')
tRNA_long$codon = str_replace_all(tRNA_long$codon,'T','U')

vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
tRNA_long$title = paste(tRNA_long$anticodon,sep="")
tRNA_long$codon = factor(tRNA_long$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
tRNA_long$title = factor(tRNA_long$title,levels= tapply(tRNA_long$title, as.integer(tRNA_long$codon),unique))

tRNA_long$color = sapply(tRNA_long$codon,function(x)substr(x,3,3))
tRNA_long$color = paste("NN",tRNA_long$color,sep="")
set_color = c(NNA="#B2DF8A",NNU="#33A02C",NNC="#1F78B4",NNG="#A6CEE3")
# %>% filter(count == 0)  # #NOT WORKNG, don't know why, 
# will use a ifelse in the size argument (beurk)
# found a better solution with stroke
# Plot it


b <- ggplot(tRNA_long[!tRNA_long$aa_name %in% c("Ter"),], aes(x=title, y=sp,color=color)) + 
  #geom_point(aes(size= ifelse(count==0, NA, count))) +  #not good to introduce NAs
  geom_point(aes(size=count), stroke = 0) + # stroke=0 so that we don't see the 0, but they are still here 
  #scale_size(range=c(0, 10)) +
  scale_size_area("Gene copy number") +                       # better eye rendering than radius
  theme_bw() +
  theme(
    axis.text.x = element_text(size=5, angle=90),
    axis.text.y=element_text(size=3),
    panel.grid.major = element_line(color = "grey",
                                    size = 0.1,
                                    linetype = 2),
    panel.grid.major.x = element_blank() ,
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "grey", fill = NA)) + theme(
      title =  element_text(size=20, family="economica"),
      legend.text =  element_text(size=20, family="economica"),
      strip.text = element_text(size=15, family="economica",face="bold"),
      axis.text.x =  element_text( size=5,vjust = 0.5, family="economica"),
      axis.text.y =  element_text( size=3, family="economica",face="italic"),
      plot.title = element_text(hjust = 0.5,margin = margin(0,0,20,0))
    ) +
  xlab(NULL) + ylab(NULL)+scale_color_manual("Decoded codons",values = set_color) + facet_wrap(~amino_acid,ncol=length(unique(tRNA_long$amino_acid)),scales="free_x")+ 
  guides(colour = guide_legend(override.aes = list(size=5)))
b

# ggMarginal(b, type="histogram")    DOES NOT WORK ;(

#p + geom_facet(panel='tRNA',    #NOT WRKING
#               data=tRNA_long, 
#               geom=geom_point,
#               mapping=aes(x=tRNA, y=sp, size=count))

# histogram of the total number of tRNA loci per species
h <- tRNA_long %>% group_by(sp) %>%
  summarize(tRNAlociCount = sum(count),clade_group=unique(clade_group)) %>% 
  ggplot(aes(y=sp, x=tRNAlociCount,fill=clade_group)) +
  geom_col(col="black",size=0.1) +
  scale_fill_manual("Clades",values=Clade_color[unique(edge_clade)]) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(size=7),
    panel.grid.major = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) + theme(
      title =  element_text(size=20, family="economica"),
      legend.text =  element_text(size=15, family="economica"),
      strip.text = element_text(size=15, family="economica",face="bold"),
      axis.text.x =  element_text( size=10, family="economica"),
      axis.text.y =  element_text( size=0, family="economica"),
      plot.title = element_text(hjust = 0.5,margin = margin(0,0,20,0))
    )+
  xlab(NULL) + ylab(NULL)
h
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#fdfd99" "#e2cc1a"
#make hist

pB = b %>% insert_left(pA, width=0.7) %>% insert_right(h, width=0.4 )
pB

jpeg(paste(path_pannel,"p2A_supp.jpg",sep=""), width = 2000/.2,  1200/.2,res=100/.2)
print(pB)
dev.off()


# Supplementary Figure 2

imgA = load.image(paste(path_pannel,"p2A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure2_supp.pdf",sep=""), width=6.75, height=4)

  m=matrix(rep(1,10*10), nrow=10)

  layout(m)

  par(mar=c(0, 1, 0, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}

