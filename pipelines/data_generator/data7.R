# Generate Data 7
library(ggplot2)
library(stringr)
library(stringi)

# Define a function to check if a string contains only A, T, C, or G
check_string <- function(string) {
  allowed_letters <- c("A", "T", "C", "G")
  all_chars <- unlist(strsplit(string, ""))
  all(all_chars %in% allowed_letters )
}

path = "/home/fbenitiere/data/"

clade_dt = read.delim(paste( "data/clade_dt.tab",sep=""),header=T)
rownames(clade_dt) = clade_dt$species


data7 = data.frame()
for (species in clade_dt$species){print(species)
  if (file.exists(paste(path , "Projet-SplicedVariants/Annotations/",species,"/GC_content.tab",sep=""))){
    dt = read.delim(paste(path , "Projet-SplicedVariants/Annotations/",species,"/GC_content.tab",sep=""))
    dt[is.na(dt$genome_character),"genome_character"] = "NA"
    
    dt = dt[sapply(dt$genome_character,check_string),] # Enleve N
    
    rownames(dt) = dt$genome_character
    total_count = tapply(dt$Freq , nchar(dt$genome_character),sum)
    dt$proba_observed = dt$Freq / total_count[nchar(dt$genome_character)]
    
    dt$proba_expected = apply(dt,1,function(x){
      dt[substr(x["genome_character"],1,1),]$proba_observed * dt[substr(x["genome_character"],2,2),]$proba_observed
    })
    
    dt$shift = dt$proba_observed - dt$proba_expected
    dt$species = species
    data7 = rbind(data7,dt)
    
  } 
}

write.table(data7,"data/data7.tab",quote=F,row.names = F,sep="\t")
