# Generate Supplementary Texts Illustration 2
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


pA = ggplot(dt,aes(x=S,y=FOP,col=as.character(lambda)))  + geom_hline(aes(yintercept=FOP0,col=as.character(lambda)),size=1.1,linetype='dashed') +
  geom_hline(aes(yintercept=1),linetype='dashed') + geom_line(size=2)+
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
  )+ scale_color_manual(substitute(paste(lambda)),values=c("#1F78B4","#33A02C")) +
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,10000),labels=c()) + 
  annotation_logticks(sides="b") + xlab("S (log scale)") + geom_text(data=data.frame(a=1),label="Optimum",x=-1,y=0.96,color="black",
                                                                     fontface= "italic", size=7, family="ubuntu condensed")+
  geom_text(data=dt[!duplicated(dt$lambda),],aes(label="No selection",y=FOP0,col=as.character(lambda)),vjust=-.5 ,x=-2.8,
            fontface= "italic", size=7, family="ubuntu condensed")+ theme(legend.position="none") + xlab("")
pA

jpeg(paste(path_pannel,"i2A_supptext.jpg",sep=""), width = 6000/4, height = 3000/4,res=600/4)
print(pA)
dev.off()

# Illustration B

pB = ggplot(dt,aes(x=S,y=logit(FOP),col=as.character(lambda)))  +
  geom_hline(aes(yintercept=logit(FOP0),col=as.character(lambda)),size=1.1,linetype='dashed') +geom_line(size=2)+
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
  )+ scale_color_manual(substitute(paste(lambda)),values=c("#1F78B4","#33A02C")) +
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,10000),labels=c(0.001,0.01,0.1,1,10,100,1000,10000)) + 
  annotation_logticks(sides="b") + xlab("S (log scale)") 
pB

jpeg(paste(path_pannel,"i2B_supptext.jpg",sep=""), width = 6000/4, height = 3000/4,res=600/4)
print(pB)
dev.off()


# Supplementary Texts Illustration 2

imgA = load.image(paste(path_pannel,"i2A_supptext.jpg",sep="") )
imgB = load.image(paste(path_pannel,"i2B_supptext.jpg",sep="") )

{
  pdf(file= paste(path_figure,"Illustration2_supptext.pdf",sep=""), width=4.5, height=4)
  
  m = matrix(rep(NA,1*2), nrow=2)
  
  m[,1]=c(1,2)
  layout(m)
  
  par(mar=c(0, 2, 1, 3.2))
  plot(imgA, axes=FALSE)
  par(mar=c(0.2, 0, 0, 0))
  plot(imgB, axes=FALSE)
  dev.off()
}
