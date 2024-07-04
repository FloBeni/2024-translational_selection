# Generate Supplementary Figure 10
source("figure/figure_supp_generator/library_path.R")
resolution = 4

# Pannel A

data5 = read.delim("data/data5_supp.tab",comment.char = "#")
data5$gene_set = str_replace_all(data5$gene_set,"all,","genes,")
data5[data5$categorie == "POC-matching triplets (POCMT)",]$categorie = "control"
data5[data5$categorie == "Putative optimal codons (POC)",]$categorie = ""

# dt_graph = data5[ data5$species "Musca_domestica" & data5$set != "POCs",]
dt_graph = data5[  data5$set != "POCs" & data5$species %in%  c( "Bactrocera_oleae","Ceratitis_capitata","Bombus_terrestris","Ignelater_luminosus"),]
dt_graph$species = paste(str_replace_all(dt_graph$species,"_"," "),sep="")


pA = ggplot(dt_graph ,
            aes(x=fpkm ,y=freq*100,fill=paste(set,categorie),col=paste(set,categorie)))  + geom_point(alpha=0)+
  geom_line(data=dt_graph,size=2,aes(linetype=paste(set,categorie))) +
  # geom_point(data=dt_graph,pch=21,col="black",size=2)+
  scale_fill_manual(values=set_color[c(2,1,4,3)]) +
  scale_color_manual(values=set_color[c(2,1,4,3)]) +
  scale_shape_manual(values=c(21,22,24,23,25,20))+
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab("Gene expression level (FPKM, log scale)") + ylab("Codon frequency (%)") + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=30,family="ubuntu condensed",vjust=0),
    axis.title.y = element_text(color="black", size=30, family="ubuntu condensed",vjust=1.7),
    axis.text.y =  element_text(color="black", size=25, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=25, family="ubuntu condensed"),
    title =  element_text(color="black", size=25, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=25, family="ubuntu condensed",margin = margin(t = 5)),
    strip.text = element_text(size=15,face="italic")
  )+labs(fill='Categories',color='Categories',shape='',linetype='')+
  guides(fill = guide_legend(override.aes = list(pch=NA),order = 1),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2),
  ) + scale_x_log10(
    breaks=c(0.01,0.1,1,10,100,1000,10000,50000),
    labels=c(0.01,0.1,1,10,100,1000,10000,50000),limits=c(0.005,1000))+ ylim(0.2*100,0.65*100) +
  annotation_logticks(sides = "b")+  
  guides(fill="none",linetype="none",shape="none") +
  facet_wrap(~species)
pA

jpeg(paste(path_pannel,"p10A_supp.jpg",sep=""),  width = 10500/2*0.8/2,  5400/2/2,res=600/1.8/2)
print(pA)
dev.off()



# Supplementary Figure 10

imgA = load.image(paste(path_pannel,"p10A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure10_supp.pdf",sep=""), width=9*0.8, height=5)
  
  m = matrix(rep(NA,100*100), nrow=100)
  
  m
  par(mar=c(1, 1, 2, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}

