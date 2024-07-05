# Generate Supplementary Figure 9
source("figure/figure_supp_generator/library_path.R")


data1 = read.delim("data/data1_supp.tab",comment.char = "#")
rownames(data1) = data1$species
data1$clade_group = GTDrift_list_species[data1$species,]$clade_group

data1 = data1[ data1$nb_codon_not_decoded == 0  & data1$pval_aa_fpkm < 0.05 & data1$nb_genes_filtered >= 5000 ,]
data1 = data1[data1$clade_group %in% c("Diptera","Lepidoptera"),]
data1 = data1[data1$species != "Eumeta_japonica",]
data1$species = paste(str_replace_all(data1$species,"_"," "),sep="")


data12 = read.delim("data/data12_supp.tab",comment.char = "#")
data12$species = paste(str_replace_all(data12$species,"_"," "),", GCi=",round(data1[data12$species,]$gci,2),sep="")
data12$third_nucl = substr(data12$codon,3,3)
data12$second_nucl = substr(data12$codon,1,2)

# Pannel A

data12$species = factor(data12$species,levels=paste(data1[order(data1$gci),]$species,", GCi=",round(data1[order(data1$gci),]$gci,2),sep=""))
dt_graph = data12[data12$aa_name == "Val",]

count_dinucl = table(dt_graph$second_nucl)
count_dinucl = count_dinucl[order(count_dinucl,decreasing = T)]
dt_graph[dt_graph$second_nucl == names(count_dinucl[1]),]$third_nucl = paste(dt_graph[dt_graph$second_nucl == names(count_dinucl[1]),]$third_nucl,"1",sep="")


vector = unique(paste(dt_graph$codon))
names( vector ) = unique(paste(dt_graph$third_nucl))

names(set_color) = c("G1","C1","A1","T1","G","C","A","T")

dt_graph$third_nucl = factor(dt_graph$third_nucl,levels = c("G1","A1","T1","C1"))


pA = ggplot(dt_graph,aes(x=fpkm,y=rscu,fill=third_nucl,col=third_nucl)) +geom_line(size=1.5) +
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=30,family="ubuntu condensed",vjust=0),
    axis.title.y = element_text(color="black", size=30, family="ubuntu condensed",vjust=1.7),
    axis.text.y =  element_text(color="black", size=25, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=25, family="ubuntu condensed"),
    title =  element_text(color="black", size=25, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=25, family="ubuntu condensed",margin = margin(t = 5)),
    strip.text = element_text(size=15,face="italic")
  )+ xlab("Gene expression level (FPKM, log scale)") +
  scale_x_log10(
    breaks=c(0.01,0.1,1,10,100,1000,10000,50000),
    labels=c(0.01,0.1,1,10,100,1000,10000,50000),limits=c(0.005,1000))+ 
  scale_alpha_manual("Codons",values=c("CDS"=1,"intronic control"=.5)) +
  scale_fill_manual("Codons",values=set_color,label=vector)  +
  scale_color_manual("Valine\nsynonymous\ncodons",values=set_color,label=vector) +
  scale_shape_manual("Codons",values=c("intronic control"=24,"CDS"=21)) + 
 annotation_logticks(sides = "b")+   facet_wrap(~species) + ylab("RSCU")+labs(fill='Codons',color='Codons',shape='',linetype='')+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 1,override.aes = list(lwd=3)),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  )

pA

jpeg(paste(path_pannel,"p9A_supp.jpg",sep=""),  width = 10500/2*0.8/2,  5400/2/2,res=600/1.8/2)
print(pA)
dev.off()


# Supplementary Figure 9

imgA = load.image(paste(path_pannel,"p9A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure9_supp.pdf",sep=""), width=9*0.8, height=5)
  
  m = matrix(rep(NA,100*100), nrow=100)
  
  m
  par(mar=c(1, 1, 2, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}

