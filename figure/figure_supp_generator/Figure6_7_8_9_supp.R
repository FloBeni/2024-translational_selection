# Generate Supplementary Figure 6, 7, 8 & 9
source("figure/figure_supp_generator/library_path.R")


library(RColorBrewer)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))
set_color_class = set_color
names(set_color_class) = c("wobble -> watson-crick","watson-crick -> wobble",
                           "non-optimal -> watson-crick",
                           "non-optimal -> wobble","watson-crick -> non-optimal",
                           "wobble -> non-optimal",
                           "optimal -> non-optimal","non-optimal -> optimal")
set_color_class[c(1,7,8)] = c("#e2cc1a","#E31A1C","#33A02C")
set_shape_class = c(CDS=21,intron=24,CDS_simulation=23)
set_linetype_class = c(CDS="solid",intron="dashed",CDS_simulation="dotted")
set_alpha_class = c(CDS=1,intron=0.5)

for (region in c("CDS","Intron")){
  for (ylabel in c( "SNP density","Subsitutions rate")){
    method_to_calculate = "per_gene"
    
    # Supplementary Pannel B
    
    if (ylabel == "SNP density"){
      data_variant = read.delim("data/data5_supp.tab") 
    }else if (ylabel == "Subsitutions rate"){ 
      data_variant = read.delim("data/data4_supp.tab")}
    
    dt_graph = data_variant[ data_variant$method_to_calculate == method_to_calculate & data_variant$group == "Wb_WC_notambiguous" ,]
    
    pB = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("Wb_WC_notambiguous")
    if (region == "CDS"){
      pB = pB + geom_point(aes(y=density_optimal_to_nonoptimal_codon,fill="optimal -> non-optimal",pch='CDS',alpha='CDS'),size=3)+
        geom_line(aes(y=simulation_density_optimal_to_nonoptimal_codon,color="optimal -> non-optimal",linetype='CDS_simulation',alpha='CDS'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_codon,
                          ymax=confint_high_density_optimal_to_nonoptimal_codon,color="optimal -> non-optimal"),
                      width=0.03,show.legend=FALSE) +
        
        geom_point(aes(y=density_nonoptimal_to_optimal_codon,fill="non-optimal -> optimal",pch='CDS',alpha='CDS'),size=3)+
        geom_line(aes(y=simulation_density_nonoptimal_to_optimal_codon,color="non-optimal -> optimal",linetype='CDS_simulation',alpha='CDS'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_codon,
                          ymax=confint_high_density_nonoptimal_to_optimal_codon,color="non-optimal -> optimal"),
                      width=0.03,show.legend=FALSE)
    } else {
      pB = pB + geom_point(aes(y=density_optimal_to_nonoptimal_intron,fill="optimal -> non-optimal",linetype='intron',pch='intron',alpha='intron'),size=3)+
        geom_line(aes(y=simulation_density_optimal_to_nonoptimal_intron,color="optimal -> non-optimal",linetype='CDS_simulation',alpha='intron'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_intron,
                          ymax=confint_high_density_optimal_to_nonoptimal_intron,color="optimal -> non-optimal"),
                      width=0.03,show.legend=FALSE) +
        
        geom_point(aes(y=density_nonoptimal_to_optimal_intron,fill="non-optimal -> optimal",linetype='intron',pch='intron',alpha='intron'),size=3)+
        geom_line(aes(y=simulation_density_nonoptimal_to_optimal_intron,color="non-optimal -> optimal",linetype='CDS_simulation',alpha='intron'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_intron,
                          ymax=confint_high_density_nonoptimal_to_optimal_intron,color="non-optimal -> optimal"),
                      width=0.03,show.legend=FALSE)
    }
    
    pB = pB + theme_bw() + theme(
      axis.title.x = element_text(color="black", size=35,family="economica"),
      axis.title.y = element_text(color="black", size=35, family="economica"),
      axis.text.y =  element_text(color="black", size=30, family="economica"),
      axis.text.x =  element_text(color="black", size=30, family="economica"),
      title =  element_text(color="black", size=20, family="economica"),
      legend.text =  element_text(color="black", size=20, family="economica")
    )+ xlab("Gene expression level (FPKM, log scale)") +
      scale_x_log10(limits = c(0.01,1000),
                    breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                    labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
      ylab(paste(region,ylabel)) + scale_fill_manual("",values=set_color_class) + 
      scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
      scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
      guides(fill= guide_legend(override.aes = list(pch=21,byrow = TRUE),order = 2)
      ) +   theme(legend.spacing.y = unit(.2, 'cm'))  + theme(legend.position='top') + annotation_logticks(sides="b")+
      guides(color = FALSE, size = FALSE , linetype=F, pch=F, alpha=F)
    pB
    
    jpeg(paste(path_pannel,"p6789B_supp.jpg",sep=""), width = 4800/1, height = 3000/1,res=400/1)
    print(pB)
    dev.off()
    
    
    # Supplementary Pannel C
    
    dt_graph = data_variant[ data_variant$method_to_calculate == method_to_calculate & data_variant$group == "duet_ambiguous" ,]
    
    pC = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("Duet not ambiguous (watson-crick -> watson-crick of the most abundant tRNA)") 
    
    if (region == "CDS"){
      pC = pC + geom_point(aes(y=density_optimal_to_nonoptimal_codon,fill="watson-crick -> wobble",pch='CDS',alpha='CDS'),size=3)+
        geom_line(aes(y=simulation_density_optimal_to_nonoptimal_codon,color="watson-crick -> wobble",linetype='CDS_simulation',alpha='CDS'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_codon,
                          ymax=confint_high_density_optimal_to_nonoptimal_codon,color="watson-crick -> wobble"),
                      width=0.03,show.legend=FALSE) +
        
        geom_point(aes(y=density_nonoptimal_to_optimal_codon,fill="wobble -> watson-crick",pch='CDS',alpha='CDS'),size=3)+
        geom_line(aes(y=simulation_density_nonoptimal_to_optimal_codon,color="wobble -> watson-crick",linetype='CDS_simulation',alpha='CDS'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_codon,
                          ymax=confint_high_density_nonoptimal_to_optimal_codon,color="wobble -> watson-crick"),
                      width=0.03,show.legend=FALSE)
    } else {
      pC = pC +  geom_point(aes(y=density_optimal_to_nonoptimal_intron,fill="watson-crick -> wobble",pch='intron',alpha='intron'),size=3)+
        geom_line(aes(y=simulation_density_optimal_to_nonoptimal_intron,color="watson-crick -> wobble",linetype='CDS_simulation',alpha='intron'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_intron,
                          ymax=confint_high_density_optimal_to_nonoptimal_intron,color="watson-crick -> wobble"),
                      width=0.03,show.legend=FALSE) +
        
        geom_point(aes(y=density_nonoptimal_to_optimal_intron,fill="wobble -> watson-crick",pch='intron',alpha='intron'),size=3)+
        geom_line(aes(y=simulation_density_nonoptimal_to_optimal_intron,color="wobble -> watson-crick",linetype='CDS_simulation',alpha='intron'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_intron,
                          ymax=confint_high_density_nonoptimal_to_optimal_intron,color="wobble -> watson-crick"),
                      width=0.03,show.legend=FALSE) 
    }
    
    pC = pC + theme_bw() + theme(
      axis.title.x = element_text(color="black", size=35,family="economica"),
      axis.title.y = element_text(color="black", size=35, family="economica"),
      axis.text.y =  element_text(color="black", size=30, family="economica"),
      axis.text.x =  element_text(color="black", size=30, family="economica"),
      title =  element_text(color="black", size=20, family="economica"),
      legend.text =  element_text(color="black", size=20, family="economica")
    )+ xlab("Gene expression level (FPKM, log scale)") +
      scale_x_log10(limits = c(0.01,1000),
                    breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                    labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
      ylab(paste(region,ylabel)) + scale_fill_manual("",values=set_color_class) + scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
      scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
      guides(fill= guide_legend(override.aes = list(pch=21,byrow = TRUE),order = 2)
      ) +   theme(legend.spacing.y = unit(.2, 'cm'))  + theme(legend.position='top') + annotation_logticks(sides="b")+
      guides(color = FALSE, size = FALSE , linetype=F, pch=F, alpha=F)
    pC
    
    jpeg(paste(path_pannel,"p6789C_supp.jpg",sep=""), width = 4800/1, height = 3000/1,res=400/1)
    print(pC)
    dev.off()
    
    
    # Supplementary Pannel D
    dt_graph = data_variant[ data_variant$method_to_calculate == method_to_calculate & data_variant$group == "ic_abondant" ,]
    
    pD = ggplot(dt_graph,aes(x=fpkm)) + ggtitle("IC (watson-crick -> watson-crick of the most abundant tRNA)") 
    
    if (region == "CDS"){
      pD = pD + geom_point(aes(y=density_optimal_to_nonoptimal_codon,fill="watson-crick -> wobble",pch='CDS',alpha='CDS'),size=3)+
        geom_line(aes(y=simulation_density_optimal_to_nonoptimal_codon,color="watson-crick -> wobble",linetype='CDS_simulation',alpha='CDS'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_codon,
                          ymax=confint_high_density_optimal_to_nonoptimal_codon,color="watson-crick -> wobble"),
                      width=0.03,show.legend=FALSE) +
        
        
        
        geom_point(aes(y=density_nonoptimal_to_optimal_codon,fill="wobble -> watson-crick",pch='CDS',alpha='CDS'),size=3)+
        geom_line(aes(y=simulation_density_nonoptimal_to_optimal_codon,color="wobble -> watson-crick",linetype='CDS_simulation',alpha='CDS'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_codon,
                          ymax=confint_high_density_nonoptimal_to_optimal_codon,color="wobble -> watson-crick"),
                      width=0.03,show.legend=FALSE) 
    } else {
      
      pD = pD + geom_point(aes(y=density_optimal_to_nonoptimal_intron,fill="watson-crick -> wobble",pch='intron',alpha='intron'),size=3)+
        geom_line(aes(y=simulation_density_optimal_to_nonoptimal_intron,color="watson-crick -> wobble",linetype='CDS_simulation',alpha='intron'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_optimal_to_nonoptimal_intron,
                          ymax=confint_high_density_optimal_to_nonoptimal_intron,color="watson-crick -> wobble"),
                      width=0.03,show.legend=FALSE) +
        
        geom_point(aes(y=density_nonoptimal_to_optimal_intron,fill="wobble -> watson-crick",pch='intron',alpha='intron'),size=3)+
        geom_line(aes(y=simulation_density_nonoptimal_to_optimal_intron,color="wobble -> watson-crick",linetype='CDS_simulation',alpha='intron'),size=1) +
        geom_errorbar(aes(ymin=confint_low_density_nonoptimal_to_optimal_intron,
                          ymax=confint_high_density_nonoptimal_to_optimal_intron,color="wobble -> watson-crick"),
                      width=0.03,show.legend=FALSE) 
    }
    
    pD = pD + theme_bw() + theme(
      axis.title.x = element_text(color="black", size=35,family="economica"),
      axis.title.y = element_text(color="black", size=35, family="economica"),
      axis.text.y =  element_text(color="black", size=30, family="economica"),
      axis.text.x =  element_text(color="black", size=30, family="economica"),
      title =  element_text(color="black", size=20, family="economica"),
      legend.text =  element_text(color="black", size=20, family="economica")
    )+ xlab("Gene expression level (FPKM, log scale)") +
      scale_x_log10(limits = c(0.01,1000),
                    breaks=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000),
                    labels=c(0.005,0.01,0.1,0.5,1,5,10,50,100,1000,10000,50000)) +
      ylab(paste(region,ylabel)) + scale_fill_manual("",values=set_color_class) + scale_color_manual("",values=set_color_class) + scale_linetype_manual("",values=set_linetype_class) +
      scale_shape_manual("",values=set_shape_class) + scale_alpha_manual("",values=set_alpha_class) +
      guides(fill= guide_legend(override.aes = list(pch=21,byrow = TRUE),order = 2)
      ) +   theme(legend.spacing.y = unit(.2, 'cm'))  + theme(legend.position='top') + annotation_logticks(sides="b")+
      guides(color = FALSE, size = FALSE , linetype=F, pch=F, alpha=F)
    pD
    
    jpeg(paste(path_pannel,"p6789D_supp.jpg",sep=""), width = 4800/1, height = 3000/1,res=400/1)
    print(pD)
    dev.off()
    
    # Figure 6, 7, 8 & 9
    
    imgB = load.image(paste(path_pannel,"p6789B_supp.jpg",sep="") )
    imgC = load.image(paste(path_pannel,"p6789C_supp.jpg",sep="") )
    imgD = load.image(paste(path_pannel,"p6789D_supp.jpg",sep="") )
    
    if (ylabel == "SNP density"){
      imgA = load.image(paste(path_require,"polymorphism.png",sep="") )
      if (region == "CDS"){
        pdf(file= paste(path_figure,"Figure6_supp.pdf",sep=""), width=10, height=7)
      } else if (region == "Intron"){
        pdf(file= paste(path_figure,"Figure8_supp.pdf",sep=""), width=10, height=7)
      }
    } else if (ylabel == "Subsitutions rate"){ 
      imgA = load.image(paste(path_require,"substitutions.png",sep="") )
      if (region == "CDS"){
        pdf(file= paste(path_figure,"Figure7_supp.pdf",sep=""), width=10, height=7)
      } else if (region == "Intron"){
        pdf(file= paste(path_figure,"Figure9_supp.pdf",sep=""), width=10, height=7)
      }
    }
    {
      
      m = matrix(rep(NA,10*10), nrow=10)
      for(i in 1:5){
        m[,i]=c(rep(1,5),rep(3,5))
      }
      for(i in 6:10){
        m[,i]=c(rep(2,5),rep(4,5))
      }
      layout(m)
      m
      
      par(mar=c(0, 0, 0, 0))
      plot(imgA, axes=FALSE)
      mtext("A",at=100,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
      par(mar=c(0, 0, 0, 0))
      plot(imgB, axes=FALSE)
      mtext("B",at=100,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
      par(mar=c(0, 0, 3, 0))
      plot(imgC, axes=FALSE)
      mtext("C",at=100,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
      par(mar=c(0, 0, 0, 0))
      plot(imgD, axes=FALSE)
      mtext("D",at=100,adj=-1, side=2, line=1, font=2, cex=1.2,las=2)
      dev.off()
    }
  }
}