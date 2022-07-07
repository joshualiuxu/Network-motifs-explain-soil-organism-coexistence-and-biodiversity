#indicator

#import triad association data and environmental data
Env= read.csv("Env.csv", row.names=1)  
Micro1= read.csv("Micro1.csv", row.names=1)  


#calculate the indicator value using TITAN2 R package
library(TITAN2)
Micro1<- Micro1[length(Micro1 >= 0)>= 3,]
Micro1<-t(Micro1)
diversity_titan<-titan(Env[,1], Micro1,numPerm=250,
boot=TRUE,nBoot=500,imax=FALSE,ivTot=FALSE,
pur.cut=0.95,rel.cut=0.95,memory=TRUE)
plot_sumz(diversity_titan,filter=TRUE,xlab=expression("Community biodiversity"),pch1=20,pch2=20, col1="#fdbb2d",col2="#1E9600")

re_diversity=diversity_titan$sppmax
write.csv(re_diversity,"re_diversity.csv")

#summary the TITAN2 result for the next forest plot
forest = read.csv("forest.csv", row.names=1)  

#forest plot
library(grid)
library(forestploter)
dt =forest
dt$'Change Point' <- paste(rep(" ", 20), collapse = " ")

dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                             sprintf("%.2f (%.2f to %.2f)",
                                     dt$est, dt$low, dt$hi))

tm <- forest_theme(base_size = 12,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

 p <- forest(dt[,c(1:5,10:12)],
            est = dt$est,
            lower = dt$low, 
            upper = dt$hi,
            ci_column = 7,
            ref_line = 4.38,
            xlim = c(3.7, 5.1),
            theme = tm)
p

ggplot2::ggsave(filename = "forest.pdf", plot = p,
                dpi = 300,
                width = 10, height = 6, units = "in")


#indicator responding to the change of biodiversity
indicator= read.csv("indicator.csv", row.names=1) 
model1 <- piecewise.linear(x = indicator$Acidobacteria, y = indicator$shannon_total, 
                          CI = TRUE, bootstrap.samples = 1000, sig.level = 0.05)
 
pdf("ind1.pdf",height=5,width=5)
plot(model1, xlab = 'The motifs from Acidobacteria', ylab = 'Community Diversity')
dev.off()

#repeat