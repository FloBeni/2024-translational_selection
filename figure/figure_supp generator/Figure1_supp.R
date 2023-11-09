source("figure/figure_supp generator/library_path.R")

############## Supplementary Pannel 1 A
data12 = read.delim("data/data12.tab")
data12$clade_group = clade_dt[data12$species,]$clade_group

data12 = data12[data12$type_aa == "Wb_WC_notambiguous",]
dt_graph = data12[ data12$nb_genes > 5000 ,]

ylabel = "var_gc3"
xlabel = "var_gci"

dt_graph[,c(ylabel,xlabel)] = sqrt(dt_graph[,c(ylabel,xlabel)])

spearman_method_aa = cor.test( dt_graph$var_gc3, dt_graph$var_gci,method="spearman",exact=F)

p1A = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=16, family="economica"),
    legend.title=element_text(size=20),
    legend.text =  element_text(color="black", size=20, family="economica",vjust = 1.5,margin = margin(t = 7)),
  ) + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Standard deviation per gene GC3") +
  xlab("Standard deviation per gene GCi") +
  ggtitle(paste(
    "rho = " ,round(spearman_method_aa$estimate,3),
    ", p-value = " ,formatC(spearman_method_aa$p.value, format = "e", digits = 1), sep="")) 
# + theme(legend.position='none')
print(p1A)

jpeg(paste(path_pannel,"p1A_supp.jpg",sep=""), width = 6000/1, height = 3500/1,res=600/1)
print(p1A)
dev.off()


############## Supplementary Pannel 1 B

ylabel = "rho_gc3_gci"
xlabel = "var_gci"

spearman_method_aa = cor.test( dt_graph$rho_gc3_gci, dt_graph$var_gci,method="spearman",exact=F)

p1B = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,fill="clade_group",label="species")) +
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black", size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=20, family="economica"),
    title =  element_text(color="black", size=16, family="economica"),
    legend.title=element_text(size=20),
    legend.text =  element_text(color="black", size=20, family="economica",vjust = 1.5,margin = margin(t = 7)),
  ) + scale_fill_manual("Clades",values=Clade_color) +
  ylab("Spearmann Rho GC3 GCi per species") +
  xlab("Standard deviation per gene GCi") +
  ggtitle(paste(
    "rho = " ,round(spearman_method_aa$estimate,3),
    ", p-value = " ,formatC(spearman_method_aa$p.value, format = "e", digits = 1), sep=""))  + theme(legend.position='none')
print(p1B)

jpeg(paste(path_pannel,"p1B_supp.jpg",sep=""), width = 4500/1, height = 3500/1,res=600/1)
print(p1B)
dev.off()




############## Supplementary Pannel 1 C
p1C = ggplot(dt_graph,aes(y=rho_gc3_gci,x=clade_group,fill=clade_group,label=species))  +
  geom_boxplot(alpha=.1) + 
  geom_hline(size=1,linetype="dashed",col="red",   yintercept = 0 ) + 
  geom_point(aes(fill=clade_group),size=3,pch=21,alpha=0.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
    axis.title.y = element_text(color="black", size=25, family="economica"),
    axis.text.y =  element_text(color="black", size=20, family="economica"),
    axis.text.x =  element_text(color="black",vjust=.5, size=0,angle = 50, family="economica"),
    title =  element_text(color="black", size=15, family="economica"),
    legend.text =  element_text(color="black", size=20, family="economica")
  ) + theme(legend.position='none') + scale_fill_manual(values=Clade_color) + 
  ylab("Spearmann Rho GC3 GCi per species")  + xlab("")  + ylim(-0.2,1) 

p1C = ggMarginal(p1C, type="histogram",fill=set_color[1]) 
p1C

jpeg(paste(path_pannel,"p1C_supp.jpg",sep=""), width = 5500/1, height = 3000/1,res=500/1)
print(p1C)
dev.off()

############## Supplementary Figure 1

imgA = load.image(paste(path_pannel,"p1A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p1B_supp.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p1C_supp.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Figure1_supp.pdf",sep=""), width=3, height=5)
  
  m = matrix(rep(c(1,2,3),1*2), nrow=3)
  
  layout(m)
  m
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  mtext("A", side=2,at=111, line=1, font=2, cex=1,las=2)
  par(mar=c(0, 0, 0, 2))
  plot(imgB, axes=FALSE)
  mtext("B", side=2,adj=-1.3,at=0, line=1, font=2, cex=1,las=2)
  par(mar=c(0, 2, 0, 0))
  plot(imgC, axes=FALSE)
  mtext("C", side=2,at=111, line=1, font=2, cex=1,las=2)
  dev.off()
}
