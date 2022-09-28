#Relationships between key factors and motifs
#ecological niche theory

#import abundance data
otu1=read.csv("otu_net1.csv",row.names=1)
otu2=read.csv("otu_net2.csv",row.names=1)
otu3=read.csv("otu_net3.csv",row.names=1)
otu4=read.csv("otu_net4.csv",row.names=1)
otu5=read.csv("otu_net5.csv",row.names=1)
otu6=read.csv("otu_net6.csv",row.names=1)
otu7=read.csv("otu_net7.csv",row.names=1)
otu8=read.csv("otu_net8.csv",row.names=1)
otu9=read.csv("otu_net9.csv",row.names=1)
otu10=read.csv("otu_net10.csv",row.names=1)
otu11=read.csv("otu_net11.csv",row.names=1)
otu12=read.csv("otu_net12.csv",row.names=1)
otu13=read.csv("otu_net13.csv",row.names=1)
otu14=read.csv("otu_net14.csv",row.names=1)
otu15=read.csv("otu_net15.csv",row.names=1)
otu16=read.csv("otu_net16.csv",row.names=1)
otu17=read.csv("otu_net17.csv",row.names=1)
otu18=read.csv("otu_net18.csv",row.names=1)
otu19=read.csv("otu_net19.csv",row.names=1)
otu20=read.csv("otu_net20.csv",row.names=1)

#save all abundance data as otulist
otulist=list(otu1,otu2,otu3,otu4,otu5,otu6,otu7,otu8,otu9,otu10,
otu11,otu12,otu13,otu14,otu15,otu16,otu17,otu18,otu19,otu20)

#rename otulist
for (i in 1:length(otulist)){
names(otulist)[i] <- paste("net", i, sep = "")
}

#calculate niche breath
#take TSEA as an example of ecological amplitude
data=c()
for(j in 1:length(otulist)){

otu=otulist[[j]]

niche=c()

for(i in 1:ncol(otu)){
envnew=env[row.names(otu),]
value = weighted.mean(envnew[,"TSEA"],otu[,i])
niche=c(niche,value)
}

niche=as.data.frame(niche)
niche$net= paste("net", j, sep = "")
data=rbind(data,niche)

print(j)
}
write.csv(data,"niche.tsea.csv")


#visualization
p1=ggplot(niche.tsea,aes(x=cycfac, y=niche))+
          geom_point(alpha=.7, size=2,color="#223e83")+
          ggpubr::stat_cor(label.y=10000)+
          ggpubr::stat_regline_equation(label.y=12000)+
          geom_smooth(aes(cycfac,niche), method=lm, se=T,color="#e71f19")+
          geom_smooth(aes(cycfac,niche), method=loess, se=T,color="#33a02c")+
          theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),panel.border = element_rect(colour = "black", fill=NA, size=1))
p1
ggsave( "tsea.cycfac.pdf",p1,width=5,height=3)


p2=ggplot(niche.tsea,aes(x=facmcom, y=niche))+
          geom_point(alpha=.7, size=2,color="#223e83")+
          ggpubr::stat_cor(label.y=10000)+
          ggpubr::stat_regline_equation(label.y=12000)+
          geom_smooth(aes(facmcom,niche), method=lm, se=T,color="#e71f19")+
          geom_smooth(aes(facmcom,niche), method=loess, se=T,color="#33a02c")+
          theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),panel.border = element_rect(colour = "black", fill=NA, size=1))
p2
ggsave( "tsea.facmcom.pdf",p2,width=5,height=3)


p3=ggplot(niche.tsea,aes(x=tranfac, y=niche))+
          geom_point(alpha=.7, size=2,color="#223e83")+
          ggpubr::stat_cor(label.y=10000)+
          ggpubr::stat_regline_equation(label.y=12000)+
          geom_smooth(aes(tranfac,niche), method=lm, se=T,color="#e71f19")+
          geom_smooth(aes(tranfac,niche), method=loess, se=T,color="#33a02c")+
          theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),panel.border = element_rect(colour = "black", fill=NA, size=1))
p3
ggsave( "tsea.tranfac.pdf",p3,width=5,height=3)

p4=ggplot(niche.tsea,aes(x=trancom, y=niche))+
          geom_point(alpha=.7, size=2,color="#223e83")+
          ggpubr::stat_cor(label.y=10000)+
          ggpubr::stat_regline_equation(label.y=12000)+
          geom_smooth(aes(trancom,niche), method=lm, se=T,color="#e71f19")+
          geom_smooth(aes(trancom,niche), method=loess, se=T,color="#33a02c")+
          theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),panel.border = element_rect(colour = "black", fill=NA, size=1))
p4
ggsave( "tsea.trancom.pdf",p4,width=5,height=3)
 
p5=ggplot(niche.tsea,aes(x=trancomfac, y=niche))+
          geom_point(alpha=.7, size=2,color="#223e83")+
          ggpubr::stat_cor(label.y=10000)+
          ggpubr::stat_regline_equation(label.y=12000)+
          geom_smooth(aes(trancomfac,niche), method=lm, se=T,color="#e71f19")+
          geom_smooth(aes(trancomfac,niche), method=loess, se=T,color="#33a02c")+
          theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),panel.border = element_rect(colour = "black", fill=NA, size=1))
p5
ggsave( "tsea.trancomfac.pdf",p5,width=5,height=3)
