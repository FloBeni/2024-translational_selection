options(scipen=999)

library(stringr)
library(ggplot2)


path = "/home/fbenitiere/data/"
# path = "/beegfs/data/fbenitiere/"

species="Drosophila_melanogaster"

data2_supp = data.frame()
for (species in c("Homo_sapiens","Drosophila_melanogaster")){print(species)
  dt = read.delim(paste(path,"Projet-NeGA/translational_selection/GC_gap_per_window_per_species/",species,".tab",sep=""))
  # da = read.delim(paste(path,"Projet-SplicedVariants/Annotations/",species,"/formatted_data/gc3_gci_per_100bp.tab.gz",sep=""))
  # da = da[da$gene_id != "gene_id",]
  # 
  # da[, c("cds_length","gene_length" ,"start","end","pos3_sites","gc3_count","posi_sites" , "gci_count")] = lapply(da[, c("cds_length","gene_length" ,"start","end","pos3_sites","gc3_count","posi_sites" , "gci_count")],as.numeric)
  # 
  # da_reduce = da[!duplicated(da$protein) ,]
  # da_reduce = da_reduce[order(da_reduce$cds_length,decreasing = T),]
  # da_reduce = da_reduce[!duplicated(da_reduce$gene_id),]
  # 
  # da = da[da$protein %in% da_reduce$protein,]
  # 
  # dt=da
  
  
  size_windows = 100
  value_windows = floor(seq(0,max(dt$start),100)/size_windows) *size_windows+size_windows/2
  names(value_windows) = as.character(seq(0,max(dt$start),100))
  
  
  dt$group = value_windows[as.character(dt$start)]
  
  size_gene = tapply(dt$pos3_sites*3 + dt$posi_sites,dt$busco_id,sum)
  dt$length_gene = size_gene[dt$busco_id]
  
  quant = quantile(size_gene,seq(0,1,1/3))
  cut = cut(size_gene,quant,include.lowest = T,include.higher=T)
  names(cut) = names(size_gene)
  table(cut)
  dt$cut = cut[dt$busco_id]
  
  
  dg = data.frame(
    species,
    per_windows=names(tapply(dt$pos3_sites,paste(dt$from,dt$group,dt$cut),sum)),
    nb_genes=tapply(dt$protein,paste(dt$from,dt$group,dt$cut),function(x) length(unique(x))),
    median_length=tapply(dt$busco_id,dt$cut,function(x) median(size_gene[unique(x)]))[tapply(as.character(dt$cut),paste(dt$from,dt$group,dt$cut),unique)],
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
  
  data2_supp = rbind(data2_supp,dg)
  
  
  dg$from = factor(dg$from,levels = c("from_5prime" , "from_3prime"))
}

write.table(data2_supp,"data/data2_supp.tab",quote=F,row.names = F,sep="\t")

