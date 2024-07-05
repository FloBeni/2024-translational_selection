# Generate Supplementary Figure 8
source("figure/figure_supp_generator/library_path.R")

# Pannel A

code = read.delim(paste("data/standard_genetic_code.tab",sep=""),comment.char = "#")
rownames(code) = code$codon

data11 = read.delim("data/data11_supp.tab",comment.char = "#")
data11[data11$nb_species_0 == "0 %" & !is.na(data11$nb_species_0) ,]$nb_species_0 = NA

data11$amino_acid = factor(data11$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa))
vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
data11$title = paste(data11$anticodon,"\n(",data11$codon,")",sep="")
data11[data11$title == "CAU\n(AUG)",]$title = "anticodon     CAU   \n(codon)      (AUG)"


data11$codon = factor(data11$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
data11$title = factor(data11$title,levels= tapply(data11$title, as.integer(data11$codon),unique))

set_color = c(A="#B2DF8A",T="#33A02C",C="#1F78B4",G="#A6CEE3")

pA = ggplot(data11,aes(x=title,y=tRNA_gene_copy,label=nb_species_0)) + geom_boxplot(aes(fill=color),outlier.shape=NA) +
  scale_fill_manual("",values = set_color) + facet_wrap(~amino_acid, ncol = 4,scales = "free")+ geom_text(family="ubuntu condensed",size=9,aes(y = y_axis_0 + 3 ),vjust=0.1) + 
  theme_bw() + theme(
    title =  element_text(size=30, family="ubuntu condensed"),
    legend.text =  element_text(size=20, family="ubuntu condensed"),
    strip.text = element_text(size=25, family="ubuntu condensed",face="bold"),
    legend.spacing.x = unit(1, 'cm'),
    legend.position="top",
    legend.box.spacing = unit(2, "cm"),
    axis.title.y = element_text(color="black",vjust=1.5),
    axis.text.x =  element_text( size=21,vjust=0.5, family="ubuntu condensed"),
    axis.text.y =  element_text( size=25, family="ubuntu condensed"),
    plot.title = element_text(hjust = 0.5,margin = margin(0,0,20,0))
  ) + theme(legend.position='none') + ylab("tRNA gene copy number") + xlab("")+
  facetted_pos_scales(
    x = list(
      NULL
    ),
    y = list(
        scale_y_continuous(limits = c(0, 55)),
        scale_y_continuous(limits = c(0, 20)),
        scale_y_continuous(limits = c(0, 30)),
        scale_y_continuous(limits = c(0, 45)),
        scale_y_continuous(limits = c(0, 30)),
        scale_y_continuous(limits = c(0, 45)),
        scale_y_continuous(limits = c(0, 30)),
        scale_y_continuous(limits = c(0, 40)),
        scale_y_continuous(limits = c(0, 35)),
        scale_y_continuous(limits = c(0, 25)),
        scale_y_continuous(limits = c(0, 55)),
        scale_y_continuous(limits = c(0, 30)),
        scale_y_continuous(limits = c(0, 30)),
        scale_y_continuous(limits = c(0, 40)),
        scale_y_continuous(limits = c(0, 40)),
        scale_y_continuous(limits = c(0, 20)),
        scale_y_continuous(limits = c(0, 30)),
        scale_y_continuous(limits = c(0, 20)),
        scale_y_continuous(limits = c(0, 20)),
        scale_y_continuous(limits = c(0, 25))
      )
    )
pA


jpeg(paste(path_pannel,"p8A_supp.jpg",sep=""), width = 2000/1,  1800/1,res=110/1)
print(pA)
dev.off()



# Supplementary Figure 8

imgA = load.image(paste(path_require,"wobble_pairing.png",sep="") )
imgB = load.image(paste(path_pannel,"p8A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure8_supp.pdf",sep=""),width=7, height=7)
  
  m=matrix(rep(NA,10*100), nrow=100)
  
  for(i in 1:10){
    m[,i]=c(rep(1,15),rep(2,85))
  }
  
  layout(m)
  m
  
  par(mar=c(0, 0, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=60,adj=-4, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(0, 0, 0.2, 0))
  plot(imgB, axes=FALSE)
  par(mar=c(0, 1, 1, 0))
  mtext("B",at=30,adj=-.5, side=2, line=1, font=2, cex=1.7,las=2)
  dev.off()
}
