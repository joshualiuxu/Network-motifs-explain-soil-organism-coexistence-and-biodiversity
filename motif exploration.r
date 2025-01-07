#explore seven network motifs

#library multiple packages
library(reshape)
library(vegan)
library(igraph)
library(bipartite)
library(car)
library(truncnorm)
library(lattice)
library(gplots)
library(ggplot2)
library(tidyr)
library(ggridges)
library(tidyverse)



#mat is the adjacency matrix obained from the network inference
#A stern reminder:
#The criteria for filtering are closely related to motif results


#The following is the definition of seven potential motifs in the R software

############################################################################
#cycle facilitation (cycfac)
############################################################################

cycfac<-function(mat){

diag(mat)=0  
mat[!lower.tri(mat, diag = TRUE)] <- 0

matp<- mat; matp[which(matp<0)]<-0; matp[which(matp>0)]<-1
nsp<-nrow(matp)
ntrip<-0

for(i in 1:nsp){
# first, subset two species facilitated by a third
nei<-sum(matp[,i])
idnei<-as.numeric(which(matp[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each facilitated pair,
# look if they are both positively associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if(matp[idnei[k],idnei[z]]==1|matp[idnei[z],idnei[k]]==1)
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}

############################################################################
# cycle competition (cyccom)
############################################################################
cyccom<-function(mat){

diag(mat)=0  
mat[!lower.tri(mat, diag = TRUE)] <- 0

matn<- mat; matn[which(matn>0)]<-0; matn[which(matn<0)]<-1
nsp<-nrow(matn)
ntrip<-0

for(i in 1:nsp){
# first, subset two species competition by a third
nei<-sum(matn[,i])
idnei<-as.numeric(which(matn[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each competition pair,
# look if they are both negatively associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if(matn[idnei[k],idnei[z]]==1|matn[idnei[z],idnei[k]]==1)
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}


############################################################################
# facilitation-mediated competition (facmcom)
############################################################################

facmcom<-function(mat){

diag(mat)=0  
mat[!lower.tri(mat, diag = TRUE)] <- 0

nsp<-nrow(mat)
mat[which(is.na(mat))]<-0
matp<- mat; matp[which(matp<0)]<-0; matp[which(matp>0)]<-1
matn<- mat; matn[which(matn>0)]<-0; matn[which(matn<0)]<-1

matt<- mat; matt[which(matt>0)]<-1; matt[which(matt<0)]<-1
#
ntrip<-0
for(i in 1:nsp){
# first, subset two species facilitated by a third
nei<-sum(matt[,i])
idnei<-as.numeric(which(matt[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each competed pair,
# look if they are both positively associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if ((matn[idnei[z], idnei[k]] + matn[ idnei[z],i] + matn[ idnei[k],i] == 2) & 
    (matp[idnei[z], idnei[k]] + matp[ idnei[z],i] + matp[ idnei[k],i] == 1))
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}


############################################################################
# competition-mediated facilitation  (commfac)
############################################################################
commfac<-function(mat){

diag(mat)=0  
mat[!lower.tri(mat, diag = TRUE)] <- 0

nsp<-nrow(mat)
mat[which(is.na(mat))]<-0
matp<- mat; matp[which(matp<0)]<-0; matp[which(matp>0)]<-1
matn<- mat; matn[which(matn>0)]<-0; matn[which(matn<0)]<-1

matt<- mat; matt[which(matt>0)]<-1; matt[which(matt<0)]<-1
#
ntrip<-0
for(i in 1:nsp){
# first, subset two species associated by a third
nei<-sum(matt[,i])
idnei<-as.numeric(which(matt[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each competed pair,
# look if they are both  associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if ((matn[idnei[z], idnei[k]] + matn[ idnei[z],i] + matn[ idnei[k],i] == 1) & 
    (matp[idnei[z], idnei[k]] + matp[ idnei[z],i] + matp[ idnei[k],i] == 2))
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}

############################################################################
# transitive facilitation  (tranfac)
############################################################################

tranfac<-function(mat){

# 计算每行和每列的 1 的数量
row_sums <- rowSums(abs(mat))
col_sums <- colSums(abs(mat))

# 按照 1 的数量对行和列进行排序（降序）
sorted_rows <- order(row_sums, decreasing = TRUE)
sorted_cols <- order(col_sums, decreasing = TRUE)

# 重新排序矩阵
sorted_mat <- mat[sorted_rows, sorted_cols]



diag(sorted_mat)=0  
sorted_mat[!lower.tri(sorted_mat, diag = TRUE)] <- 0

nsp<-nrow(sorted_mat)
sorted_mat[which(is.na(sorted_mat))]<-0
matp<- sorted_mat; matp[which(matp<0)]<-0; matp[which(matp>0)]<-1
matn<- sorted_mat; matn[which(matn>0)]<-0; matn[which(matn<0)]<-1

matt<- sorted_mat; matt[which(matt>0)]<-1; matt[which(matt<0)]<-1
#
ntrip<-0
for(i in 1:nsp){
# first, subset two species associated by a third
nei<-sum(matt[,i])
idnei<-as.numeric(which(matt[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each competed pair,
# look if they are both  associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if ((matn[idnei[z], idnei[k]] + matn[ idnei[z],i] + matn[ idnei[k],i] == 0) & 
    (matp[idnei[z], idnei[k]] + matp[ idnei[z],i] + matp[ idnei[k],i] == 2))
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}


############################################################################
# transitive competition  (trancom)
############################################################################
trancom<-function(mat){

# 计算每行和每列的 1 的数量
row_sums <- rowSums(abs(mat))
col_sums <- colSums(abs(mat))

# 按照 1 的数量对行和列进行排序（降序）
sorted_rows <- order(row_sums, decreasing = TRUE)
sorted_cols <- order(col_sums, decreasing = TRUE)

# 重新排序矩阵
sorted_mat <- mat[sorted_rows, sorted_cols]



diag(sorted_mat)=0  
sorted_mat[!lower.tri(sorted_mat, diag = TRUE)] <- 0

nsp<-nrow(sorted_mat)
sorted_mat[which(is.na(sorted_mat))]<-0
matp<- sorted_mat; matp[which(matp<0)]<-0; matp[which(matp>0)]<-1
matn<- sorted_mat; matn[which(matn>0)]<-0; matn[which(matn<0)]<-1

matt<- sorted_mat; matt[which(matt>0)]<-1; matt[which(matt<0)]<-1
#
ntrip<-0
for(i in 1:nsp){
# first, subset two species associated by a third
nei<-sum(matt[,i])
idnei<-as.numeric(which(matt[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each competed pair,
# look if they are both  associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if ((matn[idnei[z], idnei[k]] + matn[ idnei[z],i] + matn[ idnei[k],i] == 2) & 
    (matp[idnei[z], idnei[k]] + matp[ idnei[z],i] + matp[ idnei[k],i] == 0))
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}


############################################################################
# transitive competition and facilitation (trancomfac)
############################################################################
trancomfac<-function(mat){

# 计算每行和每列的 1 的数量
row_sums <- rowSums(abs(mat))
col_sums <- colSums(abs(mat))

# 按照 1 的数量对行和列进行排序（降序）
sorted_rows <- order(row_sums, decreasing = TRUE)
sorted_cols <- order(col_sums, decreasing = TRUE)

# 重新排序矩阵
sorted_mat <- mat[sorted_rows, sorted_cols]



diag(sorted_mat)=0  
sorted_mat[!lower.tri(sorted_mat, diag = TRUE)] <- 0

nsp<-nrow(sorted_mat)
sorted_mat[which(is.na(sorted_mat))]<-0
matp<- sorted_mat; matp[which(matp<0)]<-0; matp[which(matp>0)]<-1
matn<- sorted_mat; matn[which(matn>0)]<-0; matn[which(matn<0)]<-1

matt<- sorted_mat; matt[which(matt>0)]<-1; matt[which(matt<0)]<-1
#
ntrip<-0
for(i in 1:nsp){
# first, subset two species associated by a third
nei<-sum(matt[,i])
idnei<-as.numeric(which(matt[,i]==1))

if(nei>=2){
# number of search = number of possible pairs = n!/2!(n-2)!
nos<- factorial(nei)/(2*factorial(nei-2))
counter<-0

# then, for each competed pair,
# look if they are both  associated
for(k in 1:(nei-1)){for(z in (k+1):nei){
counter<-counter+1
if ((matn[idnei[z], idnei[k]] + matn[ idnei[z],i] + matn[ idnei[k],i] == 1) & 
    (matp[idnei[z], idnei[k]] + matp[ idnei[z],i] + matp[ idnei[k],i] == 1))
ntrip <- ntrip+1
}
#if(counter==nos) break
}
}
}
return(ntrip)
}


#An integrated function which summarize all types of motifs
netmotif=function(mat){
data <- c()
data = c(data, cycfac(mat))
data = c(data, cyccom(mat))
data = c(data, facmcom(mat))
data = c(data, commfac(mat))
data = c(data, tranfac(mat))
data = c(data, trancom(mat))
data = c(data, trancomfac(mat))
data=t(as.data.frame(data))
colnames(data)=c("cycfac","cyccom","facmcom","commfac","tranfac","trancom","trancomfac")
return(data)
}


#########################################################################################
#motifs within real networks and random networks
#########################################################################################

#########################################################################################
#motifs within real networks
#########################################################################################

#take the clustering region #1 as an example
#get the number of motifs in the ecological network
#import edgelist file to R space
edge1 <- read.csv("edge1.csv")
adj_mat1 <- read.csv("adj_mat1.csv", row.names=1)

motif=c()
motif=as.data.frame(netmotif(as.matrix(adj_mat1))) 


#########################################################################################
#motifs within random networks
#########################################################################################

#create random networks (rn)
rn1=list()
for (x in seq_len(999L)) {
df=matrix(data = NA,nrow = nrow(edge1),ncol = 3,byrow = FALSE,dimnames = NULL)
df[,1] = edge1[sample(nrow(edge1), nrow(edge1)),1 ]
df[,2] = edge1[sample(nrow(edge1), nrow(edge1)),2 ]
df[,3] = edge1[sample(nrow(edge1), nrow(edge1)),4 ]
rn1[[x]]=df
}

random_motif1=c()
for (x in seq_len(999L)) {
ignew=graph_from_edgelist(as.matrix(rn1[[x]][,-3]) , directed = FALSE)  
ignew=set_edge_attr(ignew, 'weight', index = E(ignew),  as.numeric(rn1[[x]][,3]) )
a=as.matrix(as_adjacency_matrix( ignew,type = "both", names = TRUE,attr="weight"))  
re=as.data.frame(netmotif(a)) 
random_motif1=rbind(random_motif1,re)
print(x)
}

write.csv(random_motif1,"random_motif1.csv")


#overrepresentation of network motifs could be calculated by z-score as described by M&M using EXCEL
