
options( stringsAsFactors = F, scipen = 999 )
library(seqinr)
library(stringr)
library(ape)
library(patchwork)
library(png)
library(ggtree)
library(caper)
library(ggExtra)
library(phylolm)
library(imager)
library(ggplot2)
library(RColorBrewer)
library(stringi)
library(forcats)
library(scales)
library(ggh4x)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))

path_figure = "figure/figure_main/"
path_pannel = "figure/pannels/"
path_require = "figure/images_library/"


wobble_type = c("T"="G-U","C"="I-C","A"="I-A","G"="U-G")

Clade_color = c(Diptera="red",Lepidoptera="#FB9A99",Coleoptera="#e2cc1a",Hymenoptera="#ba8e18","Other Insects"="#FF7F00",
                Nematoda="#B2DF8A",Mammalia="#66281A",Aves="#5b5b5b","Other Tetrapods"="#A6CEE3",Teleostei="#1F78B4","Other Metazoans"="#f5b48a","branch"="black"
)



arbrePhylo = read.tree(paste("data/GTDrift_Metazoa_phylogenetic_tree.nwk",sep=""))

life_history_traits = read.delim("data/GTDrift_life_history_traits_and_polymorphism_derived_Ne.tab")
rownames(life_history_traits) = life_history_traits$species

GTDrift_list_species = read.delim("data/GTDrift_list_species.tab")
rownames(GTDrift_list_species) = GTDrift_list_species$species


GTDrift_list_species[GTDrift_list_species$clade_group == "Other Vertebrates",]$clade_group = "Other Tetrapods"


GTDrift_list_species$clade_group = factor(GTDrift_list_species$clade_group, levels = names(Clade_color))

GTDrift_list_species$length_cm = life_history_traits[GTDrift_list_species$species,"length_cm"]
GTDrift_list_species$lifespan_days = life_history_traits[GTDrift_list_species$species,"lifespan_days"]
GTDrift_list_species$mass_kg = life_history_traits[GTDrift_list_species$species,"mass_kg"]

listNomSpecies = tapply(GTDrift_list_species$species,GTDrift_list_species$clade_group,function(x)  str_replace_all(x,"_"," "))



fitted_model <- function(x=dt_graph[,xlabel],y=dt_graph[,ylabel],label=dt_graph$species,tree = NA,display_other=T,pagels_obliged=F,lm_obliged=F){
  # Function to choose which model between PGLS, LM and Pagel's lambda is best suited to the data.
  dt_fit = data.frame()
  if ( length(tree) != 1){
    shorebird <- comparative.data(tree,data.frame(label=label,x_shorebird=x,y_shorebird=y), label, vcv=TRUE)
    fit = pgls(y_shorebird~x_shorebird,shorebird)
    summ_fit = summary(fit)
    dt_fit = rbind(dt_fit,data.frame(
      model="PGLS",
      p_val_slope = summ_fit$coefficients[2,4],
      r.squared = summ_fit$r.squared,
      adj.r.squared = summ_fit$adj.r.squared,
      aic = AIC(fit),
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    ))
    
    fit <- phylolm(y_shorebird~x_shorebird, phy = shorebird$phy, data = shorebird$data, model = "lambda")
    summ_fit = summary(fit)
    dt_fit = rbind(dt_fit,data.frame(
      model="Pagel's ",
      p_val_slope = summ_fit$coefficients[2,4],
      r.squared = summ_fit$r.squared,
      adj.r.squared = summ_fit$adj.r.squared,
      aic = AIC(fit),
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    ))
    
    if (pagels_obliged){
      
      dt_fit = rbind(dt_fit,data.frame(
        model="Pagel's ",
        p_val_slope = summ_fit$coefficients[2,4],
        r.squared = summ_fit$r.squared,
        adj.r.squared = summ_fit$adj.r.squared,
        aic = -1000000000000,
        slope = coef(fit)[2],
        intercept = coef(fit)[1]
      ))
    }
  } 
  
  fit = lm(y~x)
  summ_fit = summary(fit)
  dt_fit = rbind(dt_fit,data.frame(
    model="LM",
    p_val_slope = summ_fit$coefficients[2,4],
    r.squared = summ_fit$r.squared,
    adj.r.squared = summ_fit$adj.r.squared,
    aic = AIC(fit),
    slope = coef(fit)[2],
    intercept = coef(fit)[1]
  ))
  
  
  if (lm_obliged){
    
    dt_fit = rbind(dt_fit,data.frame(
      model="LM",
      p_val_slope = summ_fit$coefficients[2,4],
      r.squared = summ_fit$r.squared,
      adj.r.squared = summ_fit$adj.r.squared,
      aic = -10000000000001,
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    ))
  }
  
  print(dt_fit)
  
  model_sub = dt_fit[dt_fit$aic != min(dt_fit$aic),]
  dt_fit = dt_fit[dt_fit$aic == min(dt_fit$aic),]
  
  model = paste(dt_fit$model,sep="")
  R2 = paste(round(dt_fit$r.squared, 3),sep="")
  
  if (dt_fit$p_val_slope < 1e-16){pvalue = "< 1e-16"} else {
    pvalue = paste("= ",formatC(dt_fit$p_val_slope, format = "e", digits = 0),sep="")
  }
  model_non_opti = ""
  
  
  if ( length(tree) != 1 & display_other){
    AIC = paste(" AIC = ",round(dt_fit$aic),",",sep="")
    model_non_opti = paste(paste("/ ",model_sub$model,": AIC = ",round(model_sub$aic),sep=""),collapse = " ")
  }
  return(list(model=model,aic=AIC,r2=R2,pvalue=pvalue,model_non_opti=model_non_opti,slope=dt_fit$slope,intercept=dt_fit$intercept))
} 


