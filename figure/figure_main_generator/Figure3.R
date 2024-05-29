# Generate Figure 3
source("figure/figure_main_generator/library_path.R")

# Pannel 3 B

code = read.delim(paste("data/standard_genetic_code.tab",sep=""))
rownames(code) = code$codon

data12 = read.delim("data/data12_supp.tab")
data12[data12$nb_species_0 == "0 %" & !is.na(data12$nb_species_0) ,]$nb_species_0 = NA

data12$amino_acid = factor(data12$amino_acid,levels = unique(code[order(code$nb_syn,code$anticodon),]$aa))
vect_debut = c("AT","GT","AC","GC","GG","CC","TC","AG","CG","CT","TT","AA","GA","CA","TG","TA")
vect_debut = str_replace_all(vect_debut,"T","U")
data12$title = paste(data12$anticodon,"\n(",data12$codon,")",sep="")
data12[data12$title == "CAU\n(AUG)",]$title = "anticodon     CAU   \n(codon)      (AUG)"


data12$codon = factor(data12$codon,levels =  unlist(lapply(vect_debut,function(x) paste(x,c("C","U","A","G"),sep=""))) ) 
data12$title = factor(data12$title,levels= tapply(data12$title, as.integer(data12$codon),unique))

set_color = c(A="#B2DF8A",T="#33A02C",C="#1F78B4",G="#A6CEE3")

pB = ggplot(data12,aes(x=title,y=abundance,label=nb_species_0)) + geom_boxplot(aes(fill=color),outlier.shape=NA) +
  scale_fill_manual("",values = set_color) + facet_wrap(~amino_acid,scales = "free")+ geom_text(family="ubuntu condensed",size=7,aes(y = y_axis_0 + 3 ),vjust=0.1) + 
  theme_bw() + theme(
    title =  element_text(size=30, family="ubuntu condensed"),
    legend.text =  element_text(size=10, family="ubuntu condensed"),
    strip.text = element_text(size=25, family="ubuntu condensed",face="bold"),
    legend.spacing.x = unit(1, 'cm'),
    legend.position="top",
    legend.box.spacing = unit(2, "cm"),
    axis.title.y = element_text(color="black",vjust=1.5),
    axis.text.x =  element_text( size=18,vjust=0.5, family="ubuntu condensed"),
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
pB


jpeg(paste(path_pannel,"p3B.jpg",sep=""), width = 2000/1,  1200/1,res=100/1)
print(pB)
dev.off()



# Figure 3

imgA = load.image(paste(path_require,"wobble_pairing.png",sep="") )
imgB = load.image(paste(path_pannel,"p3B.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure3.pdf",sep=""), width=7, height=4.8)
  
  m=matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:10){
    m[,i]=c(rep(1,2),rep(2,8))
  }
  
  layout(m)
  m
  
  par(mar=c(0, 0, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=60,adj=-6, side=2, line=1, font=2, cex=1.4,las=2)
  
  par(mar=c(0, 0, 0.2, 0))
  plot(imgB, axes=FALSE)
  par(mar=c(0, 1, 1, 0))
  mtext("B",at=30,adj=-1, side=2, line=1, font=2, cex=1.4,las=2)
  dev.off()
}
