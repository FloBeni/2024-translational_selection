# Generate Supplementary Figure 7
source("figure/figure_supp_generator/library_path.R")
resolution = 4

# Load data from Lynch et al. 2023
# from https://www.embopress.org/doi/suppl/10.15252/embr.202357561/suppl_file/embr202357561-sup-0002-datasetev1.xlsx to which was added C.nigoni Ne
lynch_dt = read.table("data/Lynch2023_embr202357561-sup-0002-metazoa.csv",header=T,sep="\t",dec=",")
lynch_dt$species = str_replace_all(lynch_dt$Species," ","_")
lynch_dt$species = sapply(lynch_dt$species ,function(x) paste(str_split_1(x,"_")[1],str_split_1(x,"_")[2],sep="_"))
lynch_dt$genus = sapply(lynch_dt$species ,function(x) str_split_1(x,"_")[1])
rownames(lynch_dt) = lynch_dt$species
lynch_dt$mass = as.numeric(lynch_dt$Dry.Mass.at.Maturity...kg.)
Ne_genus = tapply(lynch_dt$Ne,lynch_dt$genus,mean)
mass_genus = tapply(lynch_dt$mass,lynch_dt$genus,function(x) mean(x,na.rm=T))


data1 = read.delim("data/data1_supp.tab",comment.char = "#")
data1$Ne = lynch_dt[data1$species,]$Ne
data1$Ne_estimate = "from genus"
data1[!is.na(data1$Ne),]$Ne_estimate = "from species"
data1[is.na(data1$Ne),]$Ne = Ne_genus[sapply(data1[is.na(data1$Ne),]$species ,function(x) str_split_1(x,"_")[1])]
data1[is.na(data1$Ne),]$Ne_estimate = ""

data1$Mass_Lynch = lynch_dt[data1$species,]$mass
data1$Mass_Lynch_estimate = "from genus"
data1[!is.na(data1$Mass_Lynch),]$Mass_Lynch_estimate = "from species"
data1[is.na(data1$Mass_Lynch),]$Mass_Lynch = mass_genus[sapply(data1[is.na(data1$Mass_Lynch),]$species ,function(x) str_split_1(x,"_")[1])]
data1[is.na(data1$Mass_Lynch),]$Mass_Lynch_estimate = ""

data1$clade_group = GTDrift_list_species[data1$species,]$clade_group
data1[,c("lifespan_days","length_cm","mass_kg")] = GTDrift_list_species[data1$species,c("lifespan_days","length_cm","mass_kg")]

dnds = read.delim("data/GTDrift_Metazoa_dNdS.tab",comment.char = "#")
rownames(dnds) = dnds$species
data1[,c("dNdS")] = dnds[data1$species,c("dNdS")]

# Pannel A
dt_graph = data1

ylabel = "Ne"
xlabel = "lifespan_days"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=log10(dt_graph[,xlabel]),y=log10(dt_graph[,ylabel]),
                            label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pA =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,shape="Ne_estimate"))  +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept) +
  geom_point(aes(fill=clade_group),size=4.5,alpha=.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=0, size=28,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=28, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=26, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=26, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=22, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .405,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot",
    legend.title = element_text( size=25, family="ubuntu condensed",margin = margin(t = 20)),
    legend.spacing = unit(200, "pt")
  )+ guides(fill = guide_legend(override.aes = list(size=5))) +
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", ", italic("P"), "-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  )+ scale_fill_manual("Clades",values=Clade_color) + scale_shape_manual(expression(paste(italic("N"[e])," estimates")),values=c(24,21)) +
  ylab(expression(paste(italic("N"[e]^pi*''^mu)))) + 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000),labels=c(0.05,0.1,0.5,1,5,10,100,1000,10000)) +
  xlab("Longevity (days, log scale)")+ 
  guides(fill = guide_legend(byrow = TRUE,override.aes = list(size=5,pch=21)),
         shape = guide_legend(byrow = TRUE,override.aes = list(fill="black")))  +
  theme(legend.spacing.y = unit(-.1, 'cm'))  +
  scale_y_log10(labels=label_log(digits = 2))  +annotation_logticks(sides="lb")
pA

jpeg(paste(path_pannel,"p7A_supp.jpg",sep=""), width = 8700/resolution, height = 6000/resolution,res=900/resolution)
print(pA)
dev.off()


# Pannel B

dt_graph = data1

ylabel = "Ne"
xlabel = "length_cm"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=log10(dt_graph[,xlabel]),y=log10(dt_graph[,ylabel]),label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pB = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,shape="Ne_estimate"))    +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=4,alpha=.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=-.5, size=25,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=25, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=23, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=23, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=24, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.7,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", ", italic("P"), "-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) + scale_fill_manual(values=Clade_color) + scale_shape_manual(expression(paste(italic("N"[e])," estimates")),values=c(24,21)) +
  ylab(expression(paste(italic("N"[e]^pi*''^mu)))) + 
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000,5000),labels=c(0.01,0.1,1,10,100,1000,5000)) + xlab("Body length (cm, log scale)")+
  annotation_logticks(sides="lb")  + scale_y_log10() + guides(fill = guide_legend(override.aes = list(size=5,pch=21)),
                                                              shape = guide_legend(override.aes = list(fill="black"))) +
  scale_y_log10(labels=label_log(digits = 2)) 
pB

jpeg(paste(path_pannel,"p7B_supp.jpg",sep=""),width = 4200/2/resolution, height = 4000/2/resolution,res=650/2/resolution)
print(pB)
dev.off()


# Pannel C

dt_graph = data1

ylabel = "Ne"
xlabel = "Mass_Lynch"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 


model_to_use = fitted_model(x=log10(dt_graph[,xlabel]),y=log10(dt_graph[,ylabel]),label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pC =  ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,shape="Mass_Lynch_estimate"))   +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=4,alpha=.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=-.5, size=26,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=26, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=24, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=24, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=24, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.7,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) +
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", ", italic("P"), "-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) + scale_fill_manual("Clades",values=Clade_color) + scale_shape_manual(expression(paste(italic("N"[e])," estimates")),values=c(24,21)) +
  ylab(expression(paste(italic("N"[e]^pi*''^mu)))) +   theme(legend.position="none")+ xlab("Body mass (kg, log scale)")+
  scale_x_log10(breaks=c(10^-12,10^-8,10^-4,10^0,10^4,10^6),labels=label_log(digits = 2),limits = c(0.00000000001,10000)) + 
  annotation_logticks(sides="lb")  + scale_y_log10() + guides(fill = guide_legend(override.aes = list(size=5,pch=21)),
                                                              shape = guide_legend(override.aes = list(fill="black"))) +
  scale_y_log10(labels=label_log(digits = 2)) 
pC

jpeg(paste(path_pannel,"p7C_supp.jpg",sep=""),width = 4200/2/resolution, height = 4000/2/resolution,res=650/2/resolution)
print(pC)
dev.off()


# Pannel D

dt_graph = data1

ylabel = "Ne"
xlabel = "dNdS"
dt_graph = dt_graph[!is.na(dt_graph[,xlabel]) & !is.na(dt_graph[,ylabel]) & dt_graph$species %in% arbrePhylo$tip.label,] 

model_to_use = fitted_model(x=log10(dt_graph[,xlabel]),y=log10(dt_graph[,ylabel]),label=dt_graph$species,tree=arbrePhylo,display_other=F,pagels_obliged=T)

pD = ggplot(dt_graph,aes_string(y=ylabel,x=xlabel,shape="Ne_estimate"))   +
  geom_abline(lwd=1,slope = model_to_use$slope, intercept = model_to_use$intercept)+
  geom_point(aes(fill=clade_group),size=4,alpha=.7) + theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=-.5, size=26,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=26, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=24, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=24, family="ubuntu condensed"),
    title =  element_text(color="black", size=20, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=24, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.7,vjust=-1, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position="none")+
  labs(
    caption = substitute(paste(model,lambda," :",aic," R"^2,"= ",r2,", ", italic("P"), "-value ",pvalue,model_non_opti), model_to_use),
    title = paste("N = ",nrow(dt_graph)," species",sep="")
  ) + scale_fill_manual(values=Clade_color) + scale_shape_manual(expression(paste(italic("N"[e])," estimates")),values=c(24,21)) +
  ylab(expression(paste(italic("N"[e]^pi*''^mu)))) + xlab(expression(italic(d)[N]/italic(d)[S] ~ "(log scale)")) +
  scale_x_log10()+ annotation_logticks(sides="lb")  + 
  scale_y_log10() + guides(fill = guide_legend(override.aes = list(size=5,pch=21)),
                           shape = guide_legend(override.aes = list(fill="black"))) +
  scale_y_log10(labels=label_log(digits = 2)) 
pD

jpeg(paste(path_pannel,"p7D_supp.jpg",sep=""),width = 4200/2/resolution, height = 4000/2/resolution,res=650/2/resolution)
print(pD)
dev.off()


# Supplementary Figure 7

imgA = load.image(paste(path_pannel,"p7A_supp.jpg",sep="") )
imgB = load.image(paste(path_pannel,"p7B_supp.jpg",sep="") )
imgC = load.image(paste(path_pannel,"p7C_supp.jpg",sep="") )
imgD = load.image(paste(path_pannel,"p7D_supp.jpg",sep="") )

{
  pdf(file=paste(path_figure,"Figure7_supp.pdf",sep=""), width=10, height=8)
  
  m = matrix(rep(NA,2*100), nrow=2)
  
  m[1,]=c(rep(1,60),rep(2,40))
  m[2,]=c(rep(3,50),rep(4,50))
  layout(m)
  
  par(mar=c(0,0, 2, 0))
  plot(imgA, axes=FALSE)
  mtext("A",at=0,adj=-1, side=2, line=1, font=2, cex=2,las=2)
  
  plot(imgB, axes=FALSE)
  mtext("B",at=0,adj=0, side=2, line=1, font=2, cex=2,las=2)
  
  plot(imgC, axes=FALSE)
  mtext("C",at=0,adj=-1, side=2, line=1, font=2, cex=2,las=2)
  
  plot(imgD, axes=FALSE)
  mtext("D",at=0,adj=-1, side=2, line=1, font=2, cex=2,las=2)
  dev.off()
}
