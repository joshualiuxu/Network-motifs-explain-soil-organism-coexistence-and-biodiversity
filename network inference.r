#The pipeline for the network inference
#Take the clustered region #1 as an example to record the whole pipeline
#The other 19 regions clustered by geostatistical algorithm were calculated as the same way

setwd("/home/admin123/xu")
library(CoNetinR)
library(dplyr)
library(seqtime)

#import filtered OTU abundance table;
#filtering criterion was based on M&M description
otu1=read.csv("otu_filtered1.csv",row.names=1)

#creat a list to save the bootstrap results
result_adj1=NULL

#manually set bootstrap 
for(x in 1: 999){
{
#randomly resample from OTU dataset
ind = sample_n(as.data.frame(colnames(otu1)) , 8) #The set sampling number is the minimum for all region clusters
z=as.character( ind$`colnames(otu1)`)
data_sampled= otu1[,z]

#get the correlation relationships between taxon-taxon by Spearman method using CoNetR package
result = getNetwork(mat = data_sampled, method="spearman", permutandboot= T,norm=TRUE,  bh=T, min.occ=2, keep.filtered=T, plot=F, report.full=T, verbose=F)
scores = result$scores

#construct the adjacency matrix
adjmatrix = matrix(nrow =  nrow(otu1), ncol =  nrow(otu1))
adjmatrix[lower.tri(adjmatrix)] = scores
adjmatrix = t(adjmatrix)
adjmatrix[lower.tri(adjmatrix)] = scores

#assign the value of diagonal NA to 0
for (i in 1: nrow(otu1)){
   for (j in 1: nrow(otu1)){
    if (is.na(adjmatrix[i,j])){
       adjmatrix[i,j] = 0
    }
    #else if (abs(adjmatrix[i,j] ) < 0.6){
    #adjmatrix[i,j] = 0
    }
   }
  }

#record one cycle and save in the list file. Total 999 cycles
result_adj1[[x]]= adjmatrix
print(paste(x,"/", length(seq_len(999)), "Number"))
}

#summary the 999 bootstrap data for each association
adj_mat1= matrix(data=NA, nrow = nrow(otu1), ncol = nrow(otu1), byrow = FALSE, dimnames = NULL)
for(i in 1: nrow(otu1)){
for(j in 1: nrow(otu1)){
tem=c()
for(x in 1: 999){
tem=c(tem,result_adj1[[x]][i,j]) 
}
adj_mat1[i,j]= median(tem)
}
print(i)
}
write.csv(adj_mat1,"adj_mat1.csv")


#RMT
#following by an RMT-based approach that determines the correlation cut-off threshold in an automatic fashion. 
#Random matrix theory was initially proposed in the 1960s as a procedure to identify phase transitions associated with noise in physics and material science,
#and was later adopted for studying the behaviours of many other complex systems, including gene co-expression network construction for
#predicting gene functions and molecular ecological network construction.


#construct function, revised by microeco R package
nnsd = function(sp){
			nns <- NULL
			for(j in 2:length(sp)){
				nn <- abs(sp[j] - sp[j-1])
				nns <- c(nns, nn)
			}
			nns
		}
rmt = function(cormat,lcor=0.6, hcor=0.8){
    s <- seq(0, 3, 0.1)
    pois <- exp(-s)
    geo <- 0.5 * pi * s * exp(-0.25 * pi * s^2)
    ps <- NULL  
    for(i in seq(lcor, hcor, 0.01)){
        cormat1 <- abs(cormat)
        cormat1[cormat1 < i] <- 0  
        eigen <- sort(eigen(cormat1)$value)
        ssp <- smooth.spline(eigen, control.spar = list(low = 0, high = 3)) 
        nnsd1 <- density(nnsd(ssp$y))
        nnsdpois <- density(nnsd(pois))
        chival1 <- sum((nnsd1$y - nnsdpois$y)^2/nnsdpois$y/512)
        ps <- rbind(ps, chival1)
        print(i*100)
    }
    ps <- cbind(ps, c(seq(lcor, hcor, 0.01)))
    tc <- ps[ps[,1] == min(ps[,1]), 2]
    tc
}


tc1 <- rmt(adj_mat1)
message("The optimized COR threshold: ", tc1, "...\n")

adj_mat1[abs(adj_mat1)< tc1] = 0

#save the adjacency matrix data
write.csv(adj_mat1,"adj_mat1.tc.csv")