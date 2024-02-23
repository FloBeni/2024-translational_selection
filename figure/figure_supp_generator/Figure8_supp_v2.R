# Generate Supplementary Figure 8
source("figure/figure_supp_generator/library_path.R")

dt_chrom = read.table(paste(path,"Projet-NeGA/translational_selection/daf_drosophila_melanogaster/processed/divegence_per_chrom.tab",sep="") ,header=T )

treshold_list =  c("chr2L"="down_20000000", "chr2R"="up_7000000", "chr4"="up_0",
                   "chr3L"="down_22000000","chr3R"="up_5000000","chrX"="down_20000000")
treshold_data = data.frame(treshold = sapply(treshold_list,function(x) as.numeric(str_split(x,'_')[[1]][2])))
treshold_data$chromosome = rownames(treshold_data)
treshold_data$orientation =  sapply(treshold_list,function(x) str_split(x,'_')[[1]][1])

pA = ggplot(dt_chrom,aes(x = (start+length/2)/1000000,y = nb_snp/length)) + 
  geom_vline(data=treshold_data,aes(xintercept=treshold/1000000),linetype="dashed",col="red",size=1.5) + 
  geom_line(size=0.5,aes(y = ere_mapped_pos/length  ,col = "D.erecta mapping fraction"))+
  geom_line(size=0.5,aes(y = sim_mapped_pos/length  ,col = "D.simulans mapping fraction"))+
  geom_line(size=1,aes(y = ere_ref_same/ere_mapped_pos  ,col = "D.erecta similarity"))+
  geom_line(size=1,aes(y = sim_ref_same/sim_mapped_pos  ,col = "D.simulans similarity"))+
  # geom_line(size=1,aes(y = sim_ref_same/sim_mapped_pos  ,col = "D.simulans similarity"))+
  theme_bw() + theme(
    axis.title.x = element_text(color="black", size=35,family="economica"),
    axis.title.y = element_text(color="black", size=35, family="economica"),
    axis.text.y =  element_text(color="black", size=30, family="economica"),
    axis.text.x =  element_text(color="black", size=30, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=35, family="economica"),
    strip.text = element_text(size = 35,family="economica")
  ) + scale_color_manual("",values=set_color[4:1]) + facet_wrap(~chromosome,ncol=2,scale="free_x") +ylab("")+  
  theme(legend.spacing.y = unit(.4, 'cm'))   + 
  guides(colour = guide_legend(byrow = TRUE,override.aes = list(lwd=5))) + xlab("Windows position (Mb)") + ylab("Fraction") +
theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

jpeg(paste(path_pannel,"p8A_supp.jpg",sep=""), width = 8000/1, height = 5000/1,res=400/1)
print(pA)
dev.off()


# Supplementary Figure 8

imgA = load.image(paste(path_pannel,"p8A_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure8_supp.pdf",sep=""), width=7, height=4.5)
  
  m=matrix(rep(1,10*10), nrow=10)
  
  layout(m)
  
  par(mar=c(0, 1, 0, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}
