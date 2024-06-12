# Generate Supplementary Texts Illustration 1
source("figure/figure_supptext_generator/library_path.R")


# Illustration A

logit <- function(x){log(x/(1-x))}

S =  seq(0.001,10,0.001)
lambda = 1
dt = data.frame(S ,
                FOP = lambda/(lambda - (1-exp(-S))/(1-exp(S))),
                lambda ,
                FOP0 =  lambda/(lambda +1))

lambda = 1/2
dt = rbind(dt, data.frame(S ,
                          FOP = lambda/(lambda - (1-exp(-S))/(1-exp(S))),
                          lambda ,
                          FOP0 =  lambda/(lambda +1))
)


pA = ggplot(dt,aes(x=S,y=logit( FOP ),col=as.character(lambda)))  + geom_line(size=2)+
  theme_bw() + theme(
    axis.title.x = element_text(color="black",vjust=0, size=28,family="ubuntu condensed"),
    axis.title.y = element_text(color="black",vjust=1.5, size=28, family="ubuntu condensed"),
    axis.text.y =  element_text(color="black", size=26, family="ubuntu condensed"),
    axis.text.x =  element_text(color="black", size=26, family="ubuntu condensed"),
    title =  element_text(color="black", size=25, family="ubuntu condensed"),
    text =  element_text(color="black", size=31, family="ubuntu condensed"),
    legend.text =  element_text(color="black", size=24, family="ubuntu condensed",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.35, face= "italic", size=20, family="ubuntu condensed"),
    plot.caption.position =  "plot"
  )+ scale_color_manual(substitute(paste(lambda)),values=c("#1F78B4","#33A02C")) 

jpeg(paste(path_pannel,"i1_supptext.jpg",sep=""), width = 6000/4, height = 5000/4,res=800/4)
print(pA)
dev.off()


# Supplementary Texts Pannel 1

imgA = load.image(paste(path_pannel,"i1_supptext.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Illustration1_supptext.pdf",sep=""), width=6, height=4)
  
  par(mar=c(0, 2, 0, 0))
  plot(imgA, axes=FALSE)
  dev.off()
}