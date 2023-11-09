# Generate Data 10
options(scipen=999)

library(stringr)
library(ggplot2)



path = "/home/fbenitiere/data/"
# path = "/beegfs/data/fbenitiere/"

data1 = read.delim("data/data1.tab")
data1 = data1[ data1$pval_aa_fpkm < 0.05,]
list_species = unique( data1$species )


data10 = data.frame()
for (species in list_species){print(species)
  if (file.exists(paste(path,"Projet-NeGA/translational_selection/GC_gap_per_window_per_species/",species,".tab",sep=""))){
    dt = read.delim(paste(path,"Projet-NeGA/translational_selection/GC_gap_per_window_per_species/",species,".tab",sep=""))
    
    size_windows = 100
    value_windows = floor(seq(0,max(dt$start),100)/size_windows) *size_windows+size_windows/2
    names(value_windows) = as.character(seq(0,max(dt$start),100))
    
    
    dt$group = value_windows[as.character(dt$start)]
    
    size_gene = tapply(dt$pos3_sites*3 + dt$posi_sites,dt$busco_id,sum)
    dt$length_gene = size_gene[dt$busco_id]
    
    quant = quantile(size_gene,seq(0,1,1))
    cut = cut(size_gene,quant,include.lowest = T,include.higher=T)
    names(cut) = names(size_gene)
    table(cut)
    dt$cut = cut[dt$busco_id]
    
    dg = data.frame(
      species,
      per_windows=names(tapply(dt$pos3_sites,paste(dt$from,dt$group,dt$cut),sum)),
      nb_genes=tapply(dt$protein,paste(dt$from,dt$group,dt$cut),function(x) length(unique(x))),
      pos3_sites = tapply(dt$pos3_sites,paste(dt$from,dt$group,dt$cut),sum),
      prop_cds = tapply(dt$pos3_sites,paste(dt$from,dt$group,dt$cut),sum) / tapply(dt$pos3_sites,dt$cut,sum)[tapply(as.character(dt$cut),paste(dt$from,dt$group,dt$cut),unique)],
      gc3_count = tapply(dt$gc3_count,paste(dt$from,dt$group,dt$cut),sum),
      posi_sites = tapply(dt$posi_sites,paste(dt$from,dt$group,dt$cut),sum),
      gci_count = tapply(dt$gci_count,paste(dt$from,dt$group,dt$cut),sum),
      site_data_gap = tapply(dt$site_data_gap,paste(dt$from,dt$group,dt$cut),sum),
      mean_gap = tapply(dt$mean_gap * dt$site_data_gap,paste(dt$from,dt$group,dt$cut),function(x) sum(x,na.rm = T)) /  tapply( dt$site_data_gap,paste(dt$from,dt$group,dt$cut),function(x) sum(x,na.rm = T)) 
    )
    
    
    dg$from = sapply(dg$per_windows,function(x) strsplit(x," ")[[1]][1])
    dg$group = as.numeric(sapply(dg$per_windows,function(x) strsplit(x," ")[[1]][2]))
    dg$cut = sapply(dg$per_windows,function(x) strsplit(x," ")[[1]][3])
    
    
    data10 = rbind(data10,dg)
    
    
    dg$from = factor(dg$from,levels = c("from_5prime" , "from_3prime"))
  }
}

write.table(data10,"data/data10.tab",quote=F,row.names = F,sep="\t")

